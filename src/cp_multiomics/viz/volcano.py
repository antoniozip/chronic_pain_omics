"""Volcano plot: pooled effect size vs −log10(adjusted p-value)."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import SIG_COLORS, set_publication_style

logger = logging.getLogger(__name__)


def plot_volcano(
    pooled_df: pd.DataFrame,
    modality: str,
    out_dir: Path,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
    top_n_labels: int = 15,
    effect_col: str = "yi_pooled",
    padj_col: str = "padj_pooled",
    label_col: str = "feature_id",
) -> Path:
    """Draw a volcano plot for pooled meta-analysis results.

    Args:
        pooled_df:     DataFrame from results/meta/{modality}/pooled_effects.csv.
        modality:      Used for title and filename.
        out_dir:       Directory to write the PDF.
        padj_threshold: Horizontal significance line.
        lfc_threshold:  Vertical effect-size lines (±).
        top_n_labels:  Number of significant features to label.

    Returns:
        Path to the saved PDF.
    """
    set_publication_style()

    df = pooled_df.dropna(subset=[effect_col, padj_col]).copy()
    df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

    sig_up   = (df[padj_col] < padj_threshold) & (df[effect_col] >  lfc_threshold)
    sig_down = (df[padj_col] < padj_threshold) & (df[effect_col] < -lfc_threshold)
    ns       = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    ax.scatter(df.loc[ns,       effect_col], df.loc[ns,       "neg_log10_padj"],
               s=6, alpha=0.5, color=SIG_COLORS["ns"],       rasterized=True)
    ax.scatter(df.loc[sig_up,   effect_col], df.loc[sig_up,   "neg_log10_padj"],
               s=8, alpha=0.85, color=SIG_COLORS["sig_up"],  rasterized=True)
    ax.scatter(df.loc[sig_down, effect_col], df.loc[sig_down, "neg_log10_padj"],
               s=8, alpha=0.85, color=SIG_COLORS["sig_down"], rasterized=True)

    # Reference lines
    ax.axhline(-np.log10(padj_threshold), color="black", lw=0.6, ls="--", alpha=0.6)
    ax.axvline( lfc_threshold,            color="black", lw=0.6, ls="--", alpha=0.6)
    ax.axvline(-lfc_threshold,            color="black", lw=0.6, ls="--", alpha=0.6)

    # Labels for top significant features
    sig_df = df[sig_up | sig_down]
    sig_all = sig_df.nsmallest(top_n_labels, padj_col) if not sig_df.empty else sig_df
    for _, row in sig_all.iterrows():
        ax.annotate(
            str(row[label_col])[:20],
            xy=(row[effect_col], row["neg_log10_padj"]),
            xytext=(3, 3), textcoords="offset points",
            fontsize=5.5, alpha=0.9,
        )

    n_up   = sig_up.sum()
    n_down = sig_down.sum()
    ax.set_title(
        f"{modality.capitalize()} — pooled meta-analysis\n"
        f"Up: {n_up}  |  Down: {n_down}  |  padj<{padj_threshold}",
        pad=4,
    )
    ax.set_xlabel("Pooled effect size (yi)")
    ax.set_ylabel("-log10(adjusted p-value)")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{modality}_volcano.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    logger.info("Volcano plot → %s", out_path)
    return out_path
