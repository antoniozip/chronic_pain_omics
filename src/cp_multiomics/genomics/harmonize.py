"""Harmonize a GWAS-SSF summary-statistics file to the project common schema.

Effect sizes are the integration currency of the meta-analysis (CLAUDE.md), so
this is the point where GWAS betas/SEs enter that currency. N is attached from
the study record because the harmonised file carries no sample-size column.

Most kept files report a `beta` + `standard_error` pair directly. Some
(roughly 1 in 8 of the kept chronic-pain GWAS) report only `odds_ratio` from a
logistic model, with `standard_error` left unpopulated ("NA") and the SE
recoverable only from the `ci_upper`/`ci_lower` (95%) bounds instead. Since
log(OR) sits on the same log-odds scale as beta, these are harmonized rather
than dropped: beta = log(odds_ratio), and se is derived from the CI width when
`standard_error` itself is missing.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from cp_multiomics.genomics.resolve import SumstatsStudy

logger = logging.getLogger(__name__)

# Hard-required GWAS-SSF columns. `beta`/`standard_error` are NOT here: at
# least one of {beta, odds_ratio} must be present (checked separately below),
# and standard_error may be absent/all-NA when ci_upper/ci_lower substitute.
# The variant identifier is also resolved separately (see _IDENTIFIER_CANDIDATES):
# fully-harmonised files name it `hm_rsid` rather than a bare `rsid`.
_REQUIRED_COLUMNS = [
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "effect_allele_frequency",
    "p_value",
]

# Variant identifier columns in preference order. `rsid` first keeps the 26
# already-clean studies byte-identical; `hm_rsid` (the lifted dbSNP id in
# fully-harmonised files) recovers the 31 that lack a bare `rsid`; `variant_id`
# (chr_pos_ref_alt) is the last resort when no rsID is present at all.
_IDENTIFIER_CANDIDATES = ["rsid", "hm_rsid", "variant_id"]

# Tokens that a missing identifier can serialize to across GWAS-SSF exports.
_EMPTY_ID_TOKENS = {"", "NA", "nan", "None", "<NA>"}

_RENAME_MAP = {
    "chromosome": "chromosome",
    "base_pair_location": "position",
    "effect_allele": "effect_allele",
    "other_allele": "other_allele",
    "effect_allele_frequency": "eaf",
    "p_value": "pval",
}

_OUTPUT_COLUMNS = [
    "rsid", "chromosome", "position", "effect_allele", "other_allele",
    "beta", "se", "eaf", "pval", "n", "study_id", "cohort", "phenotype", "build",
    "effect_unit",
]

# qnorm(0.975): converts a 95% CI half-width (on the log scale) to an SE.
_Z_975 = 1.959964


def _resolve_rsid(df: pd.DataFrame, name: str) -> pd.Series:
    """Return a variant-identifier Series, preferring bare rsid over hm_rsid.

    Raises ValueError only when no candidate identifier column exists at all.
    If candidate columns exist but are entirely empty, returns an all-NaN
    Series so those rows drop downstream (consistent with missing-rsid rows).
    """
    present = [c for c in _IDENTIFIER_CANDIDATES if c in df.columns]
    if not present:
        raise ValueError(
            f"{name}: missing variant identifier, need one of {_IDENTIFIER_CANDIDATES}"
        )
    # Fill gaps rather than picking one whole column: a file whose `rsid` is
    # only partly populated would otherwise lose every row that column misses,
    # visible afterwards only as a smaller row count in the log.
    resolved = pd.Series(pd.NA, index=df.index, dtype="string")
    for col in present:
        ids = df[col].astype("string").str.strip()
        ids = ids.where(~ids.isin(_EMPTY_ID_TOKENS))
        resolved = resolved.fillna(ids)
    return resolved


def harmonize_sumstats(
    gz_path: Path,
    study: SumstatsStudy,
    out_path: Path,
    build: str = "GRCh38",
) -> int:
    """Map a harmonised .h.tsv.gz to the common schema; return rows written."""
    # `n` is consumed as the METAL SAMPLESIZE weight (metal_runner.py emits
    # `WEIGHT n`), and as the tie-break for which study represents a cohort.
    # Nothing downstream filters it, and parse_sample_size() legitimately
    # returns -1 for text it cannot parse: a negative weight either makes METAL
    # abort the group or inverts that study's contribution to the pooled Z.
    # Refusing the study here makes 04b record `harmonize_failed` with a reason.
    if isinstance(study.n, bool) or not isinstance(study.n, (int, np.integer)):
        raise ValueError(
            f"{study.accession}: sample size must be a positive integer, got {study.n!r}"
        )
    if study.n <= 0:
        raise ValueError(
            f"{study.accession}: sample size must be a positive integer, got {study.n!r}"
        )

    df = pd.read_csv(gz_path, sep="\t", compression="gzip", low_memory=False)

    missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{gz_path.name}: missing GWAS-SSF columns {sorted(missing)}")

    rsid = _resolve_rsid(df, gz_path.name)

    if "beta" in df.columns and pd.to_numeric(df["beta"], errors="coerce").notna().any():
        beta = pd.to_numeric(df["beta"], errors="coerce")
        effect_unit = "beta"
    elif "odds_ratio" in df.columns:
        odds_ratio = pd.to_numeric(df["odds_ratio"], errors="coerce")
        # Non-positive odds ratios are undefined on the log-odds scale; mask
        # them to NaN before log() so they never surface as +/-inf.
        beta = np.log(odds_ratio.where(odds_ratio > 0))
        effect_unit = "log_OR"
    else:
        raise ValueError(
            f"{gz_path.name}: missing effect column, need one of ['beta', 'odds_ratio']"
        )

    if "standard_error" in df.columns:
        se = pd.to_numeric(df["standard_error"], errors="coerce")
    else:
        se = pd.Series(np.nan, index=df.index)

    if "ci_upper" in df.columns and "ci_lower" in df.columns:
        ci_upper = pd.to_numeric(df["ci_upper"], errors="coerce")
        ci_lower = pd.to_numeric(df["ci_lower"], errors="coerce")
        # GWAS-SSF CI bounds sit on the same scale as the effect column, so the
        # conversion to an SE depends on which effect column we read. Applying
        # the log form to a beta-scale CI overstates the SE by an order of
        # magnitude (0.05-0.15 would give 0.280 instead of 0.026).
        if effect_unit == "log_OR":
            # Bounds are odds ratios: move to log-odds first. Non-positive or
            # NaN bounds are undefined there and are expected on most rows
            # (only rows with a missing standard_error use this path), so the
            # log() warnings are silenced and those rows masked out.
            usable_ci = se.isna() & (ci_upper > 0) & (ci_lower > 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                se_from_ci = (np.log(ci_upper) - np.log(ci_lower)) / (2 * _Z_975)
        else:
            # Bounds are already on the beta scale. They may legitimately be
            # negative, so no positivity guard applies here.
            usable_ci = se.isna() & ci_upper.notna() & ci_lower.notna()
            se_from_ci = (ci_upper - ci_lower) / (2 * _Z_975)
        se = se.where(~usable_ci, se_from_ci)

    out = df[list(_RENAME_MAP)].rename(columns=_RENAME_MAP)
    out["rsid"] = rsid
    out["beta"] = beta
    out["se"] = se
    out["effect_unit"] = effect_unit

    # Belt-and-suspenders: no +/-inf may ever reach the output, regardless of
    # source column or upstream masking gaps.
    out["beta"] = out["beta"].replace([np.inf, -np.inf], np.nan)
    out["se"] = out["se"].replace([np.inf, -np.inf], np.nan)

    # Drop reasons are counted separately: a study that shrinks needs to be
    # attributable, and "se <= 0" is invisible to dropna() because a negative
    # or zero SE is neither NaN nor infinite. Both then enter METAL — a negative
    # SE flips the variance sign, a zero SE is an infinite weight.
    n_total = len(out)
    no_id = out["rsid"].isna() | out["rsid"].astype(str).str.strip().eq("") | \
        out["rsid"].astype(str).eq("NA")
    out = out.loc[~no_id]
    n_no_id = n_total - len(out)

    bad_se = ~(out["se"] > 0)          # NaN, negative and zero all land here
    n_bad_se = int(bad_se.sum())
    out = out.loc[~bad_se]

    n_before_beta = len(out)
    out = out.dropna(subset=["beta"])
    n_no_beta = n_before_beta - len(out)

    if n_no_id or n_bad_se or n_no_beta:
        logger.info(
            "[%s] dropped %d/%d rows: %d without an identifier, %d with a "
            "non-finite or non-positive se, %d missing beta",
            study.accession, n_total - len(out), n_total,
            n_no_id, n_bad_se, n_no_beta,
        )

    out["n"] = study.n
    out["study_id"] = study.accession
    out["cohort"] = study.cohort
    out["phenotype"] = study.phenotype
    out["build"] = build

    out = out[_OUTPUT_COLUMNS]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    n_rows = len(out)
    if n_rows == 0:
        logger.warning("[%s] harmonized 0 variants", study.accession)
    else:
        logger.info("[%s] harmonized %d variants", study.accession, n_rows)
    return n_rows
