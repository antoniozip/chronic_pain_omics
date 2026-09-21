"""Cross-modal heatmap: top significant features × modalities."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import set_publication_style

logger = logging.getLogger(__name__)

MODALITY_ORDER = ["genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"]


def plot_cross_modal_heatmap(
    pooled_by_modality: dict[str, pd.DataFrame],
    out_dir: Path,
    top_n: int = 40,
    padj_threshold: float = 0.05,
    effect_col: str = "yi_pooled",
    padj_col: str = "padj_pooled",
    label_col: str = "feature_id",
) -> Path:
    """Heatmap of pooled effect sizes for top features across all modalities.

    Features are selected as those significant (padj < threshold) in at least
    one modality; rows = features, columns = modalities.

    Args:
        pooled_by_modality: dict mapping modality name → pooled_effects DataFrame.
        out_dir:            Output directory.
        top_n:              Maximum number of features to display.

    Returns:
        Path to the saved PDF.
    """
    set_publication_style()

    # Collect significant features per modality
    sig_features: set[str] = set()
    for mod, df in pooled_by_modality.items():
        if df.empty or padj_col not in df.columns:
            continue
        sig = df.loc[df[padj_col] < padj_threshold, label_col]
        sig_features.update(sig.tolist())

    if not sig_features:
        logger.warning("No significant features found across modalities for heatmap.")
        return out_dir / "cross_modal_heatmap_empty.txt"

    # Build matrix: features × modalities
    modalities_present = [m for m in MODALITY_ORDER if m in pooled_by_modality]
    matrix_rows: dict[str, dict[str, float]] = {f: {} for f in sig_features}

    for mod in modalities_present:
        df = pooled_by_modality[mod]
        if df.empty:
            continue
        for _, row in df.iterrows():
            feat = row[label_col]
            if feat in matrix_rows:
                matrix_rows[feat][mod] = row[effect_col]

    matrix_df = pd.DataFrame(matrix_rows).T.reindex(columns=modalities_present)
    matrix_df = matrix_df.dropna(how="all")

    # Select top_n by mean |effect| across modalities
    matrix_df["_mean_abs"] = matrix_df.abs().mean(axis=1)
    matrix_df = matrix_df.nlargest(top_n, "_mean_abs").drop(columns="_mean_abs")

    if matrix_df.empty:
        logger.warning("Matrix empty after filtering.")
        return out_dir / "cross_modal_heatmap_empty.txt"

    # Cluster rows by effect pattern (hierarchical)

    data = matrix_df.fillna(0).values
    if len(data) > 1:
        order = _cluster_order(data)
        matrix_df = matrix_df.iloc[order]

    vmax = np.nanpercentile(np.abs(matrix_df.values), 95)
    fig_h = max(4.0, len(matrix_df) * 0.18)
    fig, ax = plt.subplots(figsize=(len(modalities_present) * 1.2 + 1.5, fig_h))

    im = ax.imshow(
        matrix_df.values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax, vmax=vmax,
        interpolation="none",
    )

    ax.set_xticks(range(len(modalities_present)))
    ax.set_xticklabels([m.capitalize() for m in modalities_present], rotation=30, ha="right")
    ax.set_yticks(range(len(matrix_df)))
    ax.set_yticklabels(matrix_df.index.tolist(), fontsize=6)

    plt.colorbar(im, ax=ax, label="Pooled effect size (yi)", shrink=0.6, pad=0.02)
    ax.set_title(f"Top {len(matrix_df)} features — cross-modal effect sizes", pad=6)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "cross_modal_heatmap.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    logger.info("Cross-modal heatmap → %s", out_path)
    return out_path


def _cluster_order(data: np.ndarray) -> list[int]:
    """Return row indices in hierarchical-clustering order (complete linkage)."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    try:
        dist = pdist(np.nan_to_num(data), metric="euclidean")
        Z    = linkage(dist, method="complete")
        return list(leaves_list(Z))
    except Exception:
        return list(range(len(data)))
