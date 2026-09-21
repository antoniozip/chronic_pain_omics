"""Cross-species concordance analysis.

Joins pooled human and animal effect sizes on orthologous feature pairs,
then computes directional concordance and effect-size correlation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .hcop import OrthologMapper, OrthologRecord

logger = logging.getLogger(__name__)


@dataclass
class ConcordanceSummary:
    modality: str
    n_ortholog_pairs: int
    n_concordant: int           # same direction
    pct_concordant: float
    pearson_r: float
    pearson_p: float
    spearman_r: float
    spearman_p: float


def build_concordance_table(
    human_df: pd.DataFrame,
    animal_df: pd.DataFrame,
    animal_species: str,
    mapper: OrthologMapper,
    padj_threshold: float = 0.05,
) -> pd.DataFrame:
    """Join human and animal pooled effects on ortholog pairs.

    Args:
        human_df:       pooled_effects.csv for human-only studies (feature_id = HGNC symbol).
        animal_df:      pooled_effects.csv for animal-only studies (feature_id = animal symbol).
        animal_species: e.g. "Mus musculus".
        mapper:         OrthologMapper instance.
        padj_threshold: significance threshold for flagging pairs.

    Returns:
        DataFrame with one row per ortholog pair and concordance columns.
    """
    if animal_df.empty or "feature_id" not in animal_df.columns:
        logger.warning("Empty animal DataFrame passed to build_concordance_table.")
        return pd.DataFrame()

    animal_symbols = animal_df["feature_id"].dropna().unique().tolist()
    ortholog_map = mapper.map_to_human(animal_symbols, from_species=animal_species)

    mapped_rows: list[dict] = []
    for _, row in animal_df.iterrows():
        rec: OrthologRecord | None = ortholog_map.get(row["feature_id"])
        if rec is None or not rec.human_symbol:
            continue
        mapped_rows.append(
            {
                "animal_feature_id": row["feature_id"],
                "human_feature_id":  rec.human_symbol,
                "animal_yi":         row["yi_pooled"],
                "animal_padj":       row["padj_pooled"],
                "animal_k":          row["k"],
                "ortholog_confidence": rec.confidence,
                "animal_species":    animal_species,
            }
        )

    if not mapped_rows:
        logger.warning("No ortholog pairs found for %s.", animal_species)
        return pd.DataFrame()

    animal_mapped = pd.DataFrame(mapped_rows)

    merged = animal_mapped.merge(
        human_df[["feature_id", "yi_pooled", "padj_pooled", "k"]].rename(
            columns={
                "feature_id":   "human_feature_id",
                "yi_pooled":    "human_yi",
                "padj_pooled":  "human_padj",
                "k":            "human_k",
            }
        ),
        on="human_feature_id",
        how="inner",
    )

    if merged.empty:
        logger.warning("No overlapping ortholog pairs after joining human/animal effects.")
        return merged

    merged["concordant"] = np.sign(merged["human_yi"]) == np.sign(merged["animal_yi"])
    merged["both_significant"] = (
        (merged["human_padj"] < padj_threshold) &
        (merged["animal_padj"] < padj_threshold)
    )
    merged["delta_yi"] = merged["human_yi"] - merged["animal_yi"]

    logger.info(
        "Concordance: %d ortholog pairs | %d concordant (%.1f%%)",
        len(merged),
        merged["concordant"].sum(),
        100 * merged["concordant"].mean(),
    )
    return merged


def compute_summary(merged: pd.DataFrame, modality: str) -> ConcordanceSummary:
    """Compute correlation and concordance statistics from the merged table."""
    from scipy import stats

    valid = merged.dropna(subset=["human_yi", "animal_yi"])
    if len(valid) < 3:
        return ConcordanceSummary(
            modality=modality, n_ortholog_pairs=len(merged),
            n_concordant=0, pct_concordant=0.0,
            pearson_r=float("nan"), pearson_p=float("nan"),
            spearman_r=float("nan"), spearman_p=float("nan"),
        )

    pr, pp = stats.pearsonr(valid["human_yi"], valid["animal_yi"])
    sr, sp = stats.spearmanr(valid["human_yi"], valid["animal_yi"])

    return ConcordanceSummary(
        modality=modality,
        n_ortholog_pairs=len(merged),
        n_concordant=int(merged["concordant"].sum()),
        pct_concordant=float(100 * merged["concordant"].mean()),
        pearson_r=float(pr),
        pearson_p=float(pp),
        spearman_r=float(sr),
        spearman_p=float(sp),
    )
