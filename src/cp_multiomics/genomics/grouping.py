"""Assign harmonized GWAS studies to broad clinical pain groups and pick one
representative study per independent cohort for METAL pooling.

Implements cohort_expansion_plan.md §7.7: inverse-variance meta-analysis assumes
independent samples, so at most one study per (cohort x group) enters the pool
(largest n wins). UK Biobank appears under many phenotype spellings (self-report,
ICD10, "degree bothered by ..."); these are the SAME participants and must not be
counted repeatedly. Empty-cohort studies are each presumed an independent cohort
(the Stage-1 rule), so they are never collapsed against one another.
"""

from __future__ import annotations

import logging
from functools import cache
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GROUPS_PATH = _REPO_ROOT / "conf" / "genomics" / "pain_groups.yaml"

# UKB Catalog studies and Stage-1 presumed-biobank studies are one cohort family
# (§7.7: "Treat Pan-UKBB and any Catalog UKB study as the same cohort").
_UKB_ALIASES = {"UKB", "PRESUMED_BIOBANK"}


@cache
def _load_config(path: str) -> tuple[tuple[str, ...], tuple[tuple[str, tuple[str, ...]], ...]]:
    """Return (nonpain_patterns, ((group, patterns), ...)) from the YAML config."""
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Pain-group config not found: {cfg_path}")
    with open(cfg_path) as f:
        raw = yaml.safe_load(f) or {}
    nonpain = tuple(p.lower() for p in raw.get("nonpain", []) or [])
    groups = tuple(
        (g["group"], tuple(p.lower() for p in g["patterns"]))
        for g in raw.get("groups", []) or []
    )
    if not groups:
        raise ValueError(f"No pain groups defined in {cfg_path}")
    return nonpain, groups


def assign_pain_group(phenotype: str, config_path: Path | None = None) -> str:
    """Map a phenotype string to a pain-group slug.

    Non-pain phenotypes are checked first and return "nonpain"; the ordered group
    rules are then walked (first substring match wins); no match returns
    "unmapped".
    """
    nonpain, groups = _load_config(str(config_path or DEFAULT_GROUPS_PATH))
    text = str(phenotype).lower()
    if any(pat in text for pat in nonpain):
        return "nonpain"
    for group, patterns in groups:
        if any(pat in text for pat in patterns):
            return group
    return "unmapped"


def _cohort_family(cohort: object, accession: str) -> str:
    """UKB family for biobank studies; a per-study sentinel for empty cohorts."""
    if cohort is None or (isinstance(cohort, float) and pd.isna(cohort)):
        label = ""
    else:
        label = str(cohort).strip()
    if label in _UKB_ALIASES:
        return "UKB"
    if label == "":
        return f"INDEP::{accession}"
    return label


def select_representatives(df: pd.DataFrame, config_path: Path | None = None) -> pd.DataFrame:
    """Add group / cohort_family / representative columns per §7.7.

    Within each (group, cohort_family) the max-n study is the representative
    (tie -> smallest accession); the rest are collapsed_overlap. Rows whose group
    is nonpain/unmapped are never representative.
    """
    out = df.copy()
    out["group"] = out["phenotype"].map(lambda p: assign_pain_group(p, config_path))
    out["cohort_family"] = [
        _cohort_family(c, a) for c, a in zip(out["cohort"], out["accession"])
    ]
    out["is_representative"] = False
    out["pooled_status"] = "collapsed_overlap"

    poolable = out[~out["group"].isin({"nonpain", "unmapped"})]
    for (_group, _fam), sub in poolable.groupby(["group", "cohort_family"]):
        rep_idx = sub.sort_values(["n", "accession"], ascending=[False, True]).index[0]
        out.at[rep_idx, "is_representative"] = True
        out.at[rep_idx, "pooled_status"] = "representative"

    out.loc[out["group"] == "nonpain", "pooled_status"] = "nonpain"
    out.loc[out["group"] == "unmapped", "pooled_status"] = "unmapped"
    return out


def poolable_groups(rep_df: pd.DataFrame) -> dict[str, int]:
    """Count independent representative cohorts per group (excludes collapsed)."""
    reps = rep_df[rep_df["is_representative"]]
    return reps.groupby("group")["cohort_family"].nunique().to_dict()
