"""Cross-species concordance figure for manuscript (Step 1.7 output).

Generates:
  manuscript/figures/cross_species_concordance.pdf
    Panel A: Scatter plot (mouse vs. human effect sizes), both-sig genes labelled
    Panel B: Scatter plot (rat vs. human effect sizes), both-sig genes labelled
    Panel C: Bar chart — % concordant direction by significance stratum
    Panel D: Table inset — concordant both-significant genes with key stats
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
    species_label: str,
) -> None:
    """Scatter: animal_yi (x) vs. human_yi (y), coloured by concordance."""
    nonsig  = df[~df["both_significant"]]
    bs_disc = df[df["both_significant"] & ~df["concordant"]]
    bs_con  = df[df["both_significant"] & df["concordant"]]

    # rasterized: this cloud is ~33,000 points per panel, and one vector path
    # each made the SVG rendition 17.6 MB of 66,203 <path> elements -- a file
    # no editor opens usefully, which defeats the point of shipping vector
    # sources. Only the background cloud is rasterized; the both-significant
    # points, their outlines and every label stay vector, as in
    # pipeline/13_human_rodent_volcano.py.
    ax.scatter(nonsig["animal_yi"], nonsig["human_yi"],
               s=4, alpha=ALPHA_BG, c=COL_NONSIG, linewidths=0, zorder=1,
               rasterized=True)
    if not bs_disc.empty:
        ax.scatter(bs_disc["animal_yi"], bs_disc["human_yi"],
                   s=30, alpha=ALPHA_SIG, c=COL_DISCORD, linewidths=0.4,
                   edgecolors="white", zorder=3, label=f"Both-sig, discordant (n={len(bs_disc)})")
    if not bs_con.empty:
        ax.scatter(bs_con["animal_yi"], bs_con["human_yi"],
                   s=50, alpha=0.9, c=COL_CONCORD, linewidths=0.5,
                   edgecolors="white", zorder=4, label=f"Both-sig, concordant (n={len(bs_con)})")
        # Label concordant genes
        for _, row in bs_con.iterrows():
            ax.annotate(
                row["human_feature_id"],
                (row["animal_yi"], row["human_yi"]),
                fontsize=pt("gene"), ha="left", va="bottom",
                xytext=(3, 2), textcoords="offset points",
                color=COL_CONCORD,
            )

    lim = max(
        abs(df["animal_yi"]).quantile(0.999),
        abs(df["human_yi"]).quantile(0.999),
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
    ax.set_xlabel(f"{species_label} log₂FC (meta-analysis)", fontsize=pt("axis"))
    ax.set_ylabel("Human log₂FC (meta-analysis)", fontsize=pt("axis"))
    ax.set_title(
        f"{species_label.split(' (')[0]} ↔ Human  (n={n_total:,} pairs)\n"
        f"{pct_con:.1f}% concordant direction",
        fontsize=pt("title"), pad=8,
    )
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


def concordant_gene_table(ax: plt.Axes, df_m: pd.DataFrame, df_r: pd.DataFrame) -> None:
    """Inset table of concordant both-significant genes."""
    con_m = df_m[df_m["both_significant"] & df_m["concordant"]].copy()
    con_r = df_r[df_r["both_significant"] & df_r["concordant"]].copy()

    # Merge mouse and rat
    all_genes = set(con_m["human_feature_id"]) | set(con_r["human_feature_id"])
    rows = []
    for gene in sorted(all_genes):
        m_row = con_m[con_m["human_feature_id"] == gene]
        r_row = con_r[con_r["human_feature_id"] == gene]
        rows.append({
            "Gene": gene,
            "Mouse FC\n(animal)":
                f"{m_row['animal_yi'].values[0]:+.2f}" if not m_row.empty else "—",
            "Rat FC\n(animal)":
                f"{r_row['animal_yi'].values[0]:+.2f}" if not r_row.empty else "—",
            "Human FC": "{:+.2f}".format(
                m_row["human_yi"].values[0] if not m_row.empty
                else r_row["human_yi"].values[0]),
            "In mouse": "✓" if not m_row.empty else "",
            "In rat":   "✓" if not r_row.empty else "",
        })

    if not rows:
        ax.text(0.5, 0.5, "No concordant both-significant genes", ha="center", va="center",
                fontsize=pt("title"), transform=ax.transAxes)
        ax.axis("off")
        return

    tbl_df = pd.DataFrame(rows)
    ax.axis("off")
    tbl = ax.table(
        cellText=tbl_df.values,
        colLabels=tbl_df.columns,
        cellLoc="center", loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(pt("tick"))
    tbl.scale(1, 1.4)

    # Header styling
    for j in range(len(tbl_df.columns)):
        tbl[(0, j)].set_facecolor("#2166ac")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    # Alternating row colours
    for i in range(1, len(rows) + 1):
        bg = "#e8f0fb" if i % 2 == 0 else "white"
        for j in range(len(tbl_df.columns)):
            tbl[(i, j)].set_facecolor(bg)

    ax.set_title("Concordant both-significant genes", fontsize=pt("title"), pad=16)


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
    concordance_bars(ax_c, df_m, df_r)
    concordant_gene_table(ax_d, df_m, df_r)

    fig.suptitle(
        "Cross-species transcriptomic concordance: rodent models ↔ human chronic pain",
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
    print("Concordant both-sig genes: "
          f"{sorted(set(con_m['human_feature_id']) | set(con_r['human_feature_id']))}")


if __name__ == "__main__":
    main()
