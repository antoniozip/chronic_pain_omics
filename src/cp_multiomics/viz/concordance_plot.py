"""Cross-species concordance scatter plot: human vs animal pooled effect sizes."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .style import SIG_COLORS, set_publication_style

logger = logging.getLogger(__name__)


def plot_concordance_scatter(
    concordance_df: pd.DataFrame,
    modality: str,
    animal_species: str,
    out_dir: Path,
    padj_threshold: float = 0.05,
    top_n_labels: int = 10,
) -> Path:
    """Scatter plot of human vs animal pooled effect sizes for ortholog pairs.

    Points are colored by significance in both species.
    Pearson r and concordance % annotated in the corner.

    Args:
        concordance_df: Output of concordance.build_concordance_table().
        modality:       Modality string (used in title/filename).
        animal_species: Species label for axis.
        out_dir:        Output directory.

    Returns:
        Path to saved PDF.
    """
    set_publication_style()

    if concordance_df.empty or "human_yi" not in concordance_df.columns:
        logger.warning("[%s] Empty concordance table — skipping scatter plot.", modality)
        return out_dir / f"{modality}_concordance_empty.txt"

    df = concordance_df.dropna(subset=["human_yi", "animal_yi"]).copy()
    if df.empty:
        logger.warning("[%s] Empty concordance table — skipping scatter plot.", modality)
        return out_dir / f"{modality}_concordance_empty.txt"

    both_sig = (df["both_significant"] if "both_significant" in df.columns
                else pd.Series(False, index=df.index))
    concordant = df["concordant"] if "concordant" in df.columns else pd.Series(True, index=df.index)

    colors = np.where(
        both_sig & concordant,   SIG_COLORS["both_sig_up"],
        np.where(
            both_sig & ~concordant, "#FF7F00",
            SIG_COLORS["ns"],
        ),
    )

    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    ax.scatter(df["human_yi"], df["animal_yi"],
               c=colors, s=18, alpha=0.75, linewidths=0, rasterized=True)

    # Diagonal (perfect concordance)
    lim_val = max(df[["human_yi", "animal_yi"]].abs().max().max() * 1.1, 0.5)
    ax.plot([-lim_val, lim_val], [-lim_val, lim_val],
            color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axhline(0, color="grey", lw=0.5, alpha=0.4)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.4)
    ax.set_xlim(-lim_val, lim_val)
    ax.set_ylim(-lim_val, lim_val)

    # Label top concordant pairs (both sig, same direction)
    top = df[both_sig & concordant].nlargest(top_n_labels, "human_yi")
    for _, row in top.iterrows():
        ax.annotate(
            str(row.get("human_feature_id", row.get("feature_id", "")))[:18],
            xy=(row["human_yi"], row["animal_yi"]),
            xytext=(3, 3), textcoords="offset points",
            fontsize=5.0, alpha=0.85,
        )

    # Stats annotation
    if len(df) >= 3:
        r, p = stats.pearsonr(df["human_yi"], df["animal_yi"])
        pct = 100 * concordant.mean()
        ax.text(0.05, 0.95,
                f"r = {r:.2f}  p = {p:.2e}\nConcordant: {pct:.0f}%",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=7, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

    ax.set_xlabel("Human pooled effect size")
    ax.set_ylabel(f"{animal_species} pooled effect size")
    ax.set_title(f"{modality.capitalize()} — human vs {animal_species}\n"
                 f"ortholog pairs (n={len(df)})", pad=4)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{modality}_concordance_{animal_species.split()[1].lower()}.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    logger.info("Concordance scatter → %s", out_path)
    return out_path
