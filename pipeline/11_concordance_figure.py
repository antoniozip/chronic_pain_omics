"""Cross-species concordance figure for manuscript (Step 1.7 output).

Generates:
  manuscript/figures/cross_species_concordance.pdf
    Panel A: Scatter plot (mouse vs. human effect sizes), both-sig genes labelled
    Panel B: Scatter plot (rat vs. human effect sizes), both-sig genes labelled
    Panel C: Scatter plot (mouse vs. rat effect sizes), pairs joined through
             their shared human ortholog (ortholog.concordance.rodent_pairs)
    Panel D: Bar chart — % concordant direction by significance stratum

A fourth panel listing the concordant both-significant genes was dropped on
2026-11-16: there are none, so it drew a sentence, not a plot.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cp_multiomics.ortholog.concordance import rodent_pairs

logger = logging.getLogger(__name__)

REPO_ROOT  = Path(__file__).resolve().parent.parent
CROSS_DIR  = REPO_ROOT / "results" / "cross_species" / "transcriptomics"
OUT_PDF    = REPO_ROOT / "manuscript" / "figures" / "cross_species_concordance.pdf"

# Authored point sizes mean nothing on their own: this figure is drawn 14 in
# wide and placed at the 6.69 in text width, so LaTeX shrinks every label to
# 48% of the size written here. The axis labels were authored at 8 pt and
# reached the page at 3.8. Sizes are therefore declared as the size they should
# *render* at, and pt() scales them up by the placement ratio -- the same
# arrangement 13_human_rodent_volcano.py uses, and the same targets, so
# Figures 3 and 4 carry one type scale between them.
FIG_W_IN = 14.0
PLACED_WIDTH_IN = 6.69      # 21 - 2*2 cm text width, from the geometry package

ON_PAGE_PT = {
    "gene": 6.0,            # inline annotations, the smallest thing here
    "tick": 7.0,
    "axis": 7.5,
    "title": 8.0,
    "panel_letter": 10.0,
    "suptitle": 10.0,
}


def pt(element: str) -> float:
    """Authored point size that renders at ON_PAGE_PT[element] once placed."""
    return ON_PAGE_PT[element] * FIG_W_IN / PLACED_WIDTH_IN
OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

PADJ_SIG   = 0.05
ALPHA_BG   = 0.04
ALPHA_SIG  = 0.5
COL_CONCORD   = "#2166ac"
COL_DISCORD   = "#d6604d"
COL_NONSIG    = "#aaaaaa"


def load_concordance(fname: str) -> pd.DataFrame:
    f = CROSS_DIR / fname
    if not f.exists():
        raise FileNotFoundError(f)
    return pd.read_csv(f)


def scatter_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_label: str,
    y_label: str = "Human",
    x: str = "animal_yi",
    y: str = "human_yi",
    name_col: str = "human_feature_id",
) -> None:
    """Scatter of two pooled effects, coloured by concordance.

    Defaults draw a rodent (x) against human (y); the mouse--rat panel passes
    its own columns. ``x_label``/``y_label`` are species names, e.g.
    "Mouse (Mus musculus)"; the part before " (" titles the panel.
    """
    df = df.rename(columns={x: "_x", y: "_y"})
    nonsig  = df[~df["both_significant"]]
    bs_disc = df[df["both_significant"] & ~df["concordant"]]
    bs_con  = df[df["both_significant"] & df["concordant"]]

    # rasterized: this cloud is ~33,000 points per panel, and one vector path
    # each made the SVG rendition 17.6 MB of 66,203 <path> elements -- a file
    # no editor opens usefully, which defeats the point of shipping vector
    # sources. Only the background cloud is rasterized; the both-significant
    # points, their outlines and every label stay vector, as in
    # pipeline/13_human_rodent_volcano.py.
    ax.scatter(nonsig["_x"], nonsig["_y"],
               s=4, alpha=ALPHA_BG, c=COL_NONSIG, linewidths=0, zorder=1,
               rasterized=True)
    if not bs_disc.empty:
        ax.scatter(bs_disc["_x"], bs_disc["_y"],
                   s=30, alpha=ALPHA_SIG, c=COL_DISCORD, linewidths=0.4,
                   edgecolors="white", zorder=3, label=f"Both-sig, discordant (n={len(bs_disc)})")
    if not bs_con.empty:
        ax.scatter(bs_con["_x"], bs_con["_y"],
                   s=50, alpha=0.9, c=COL_CONCORD, linewidths=0.5,
                   edgecolors="white", zorder=4, label=f"Both-sig, concordant (n={len(bs_con)})")
        # Label concordant genes
        for _, row in bs_con.iterrows():
            ax.annotate(
                row[name_col],
                (row["_x"], row["_y"]),
                fontsize=pt("gene"), ha="left", va="bottom",
                xytext=(3, 2), textcoords="offset points",
                color=COL_CONCORD,
            )

    lim = max(
        abs(df["_x"]).quantile(0.999),
        abs(df["_y"]).quantile(0.999),
        0.5,
    ) * 1.1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.plot([-lim, lim], [-lim, lim], lw=0.7, ls=":", color="grey", alpha=0.5)

    n_total = len(df)
    n_con   = df["concordant"].sum()
    pct_con = n_con / n_total * 100 if n_total > 0 else 0
    ax.set_xlabel(f"{x_label} log₂FC (meta-analysis)", fontsize=pt("axis"))
    ax.set_ylabel(f"{y_label.split(' (')[0]} log₂FC (meta-analysis)", fontsize=pt("axis"))
    ax.set_title(
        f"{x_label.split(' (')[0]} ↔ {y_label.split(' (')[0]}  (n={n_total:,} pairs)\n"
        f"{pct_con:.1f}% concordant direction",
        fontsize=pt("title"), pad=8,
    )
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=pt("tick"), loc="upper left")
    ax.tick_params(labelsize=pt("tick"))


def concordance_bars(ax: plt.Axes, df_m: pd.DataFrame, df_r: pd.DataFrame) -> None:
    """Bar chart: % concordant direction by significance stratum."""
    def pct_concordant(df: pd.DataFrame, mask) -> float:
        sub = df[mask]
        if len(sub) == 0:
            return np.nan
        return sub["concordant"].mean() * 100

    strata = ["All pairs", "Animal sig\n(padj<0.05)",
              "Human sig\n(padj<0.05)", "Both sig\n(padj<0.05)"]
    masks_m = [
        np.ones(len(df_m), dtype=bool),
        df_m["animal_padj"] < PADJ_SIG,
        df_m["human_padj"]  < PADJ_SIG,
        df_m["both_significant"],
    ]
    masks_r = [
        np.ones(len(df_r), dtype=bool),
        df_r["animal_padj"] < PADJ_SIG,
        df_r["human_padj"]  < PADJ_SIG,
        df_r["both_significant"],
    ]

    pcts_m = [pct_concordant(df_m, m) for m in masks_m]
    pcts_r = [pct_concordant(df_r, m) for m in masks_r]
    ns_m   = [int(np.sum(m)) for m in masks_m]
    ns_r   = [int(np.sum(m)) for m in masks_r]

    x    = np.arange(len(strata))
    w    = 0.35
    bars_m = ax.bar(x - w/2, pcts_m, w, color=COL_CONCORD, alpha=0.8, label="Mouse")
    bars_r = ax.bar(x + w/2, pcts_r, w, color="#66bd63", alpha=0.8, label="Rat")

    ax.axhline(50, color="grey", ls="--", lw=0.8, label="Chance (50%)")
    ax.set_xticks(x)
    ax.set_xticklabels(strata, fontsize=pt("tick"))
    ax.set_ylabel("% concordant direction", fontsize=pt("axis"))
    ax.set_ylim(0, 100)
    ax.set_title("Direction concordance by significance stratum",
                 fontsize=pt("title"), pad=8)
    ax.legend(fontsize=pt("tick"))

    # Add n= annotations
    for bar, n in zip(bars_m, ns_m):
        ax.text(bar.get_x() + bar.get_width()/2, 2, f"n={n:,}", ha="center", va="bottom",
                fontsize=pt("gene"), rotation=90, color="white")
    for bar, n in zip(bars_r, ns_r):
        ax.text(bar.get_x() + bar.get_width()/2, 2, f"n={n:,}", ha="center", va="bottom",
                fontsize=pt("gene"), rotation=90, color="white")
    ax.tick_params(labelsize=pt("tick"))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    df_m = load_concordance("concordance_musculus.csv")
    df_r = load_concordance("concordance_norvegicus.csv")
    logger.info("Mouse: %d pairs  |  Rat: %d pairs", len(df_m), len(df_r))

    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.46, wspace=0.32)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    scatter_panel(ax_a, df_m, "Mouse (Mus musculus)")
    scatter_panel(ax_b, df_r, "Rat (Rattus norvegicus)")
    pairs = rodent_pairs(df_m, df_r)
    scatter_panel(ax_c, pairs, "Mouse (Mus musculus)",
                  y_label="Rat (Rattus norvegicus)", x="mouse_yi", y="rat_yi",
                  name_col="mouse_feature_id")
    concordance_bars(ax_d, df_m, df_r)

    fig.suptitle(
        "Cross-species transcriptomic concordance in chronic pain",
        fontsize=pt("suptitle"), fontweight="bold", y=1.01,
    )

    # Panel labels
    for ax, lbl in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        ax.text(-0.10, 1.18, lbl, transform=ax.transAxes,
                fontsize=pt("panel_letter"), fontweight="bold", va="top")

    # dpi governs the rasterized cloud only. 300 is the usual journal floor for
    # raster content; 200 was below it and visible on close inspection.
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved → %s", OUT_PDF)

    # Print summary
    con_m = df_m[df_m["both_significant"] & df_m["concordant"]]
    con_r = df_r[df_r["both_significant"] & df_r["concordant"]]
    print(f"\nMouse: {len(df_m):,} pairs, {df_m['both_significant'].sum()} both-sig, "
          f"{len(con_m)} concordant both-sig")
    print(f"Rat:   {len(df_r):,} pairs, {df_r['both_significant'].sum()} both-sig, "
          f"{len(con_r)} concordant both-sig")
    print(f"Mouse-rat: {len(pairs):,} pairs, "
          f"{pairs['concordant'].mean() * 100:.1f}% concordant direction")


if __name__ == "__main__":
    main()
