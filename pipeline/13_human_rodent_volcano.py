"""Human vs. rodent volcano comparison figure (Step 1.6 follow-up).

Generates:
  manuscript/figures/human_rodent_volcano.pdf
    Panels: one volcano per human study + pooled rodent meta-analysis volcano
    Format: 3-column grid
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import offset_copy

logger = logging.getLogger(__name__)

REPO_ROOT     = Path(__file__).resolve().parent.parent
PER_STUDY_DIR = REPO_ROOT / "results" / "per_study" / "transcriptomics"
PER_STUDY_CSV = PER_STUDY_DIR / "effects_normalized_all.csv"
POOLED_CSV    = REPO_ROOT / "results" / "meta" / "transcriptomics" / "pooled_effects.csv"
OUT_PDF       = REPO_ROOT / "manuscript" / "figures" / "human_rodent_volcano.pdf"
OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

HUMAN_STUDIES = {
    "GSE126611": "Neuropathic pain\n(blood, human)",
    "GSE177034": "Low back pain\n(blood, human)",
    "GSE221921": "Nociplastic/FM\n(blood, human)",
    "GSE250152": "Morton's neuroma\n(nerve, human)",
    "GSE235287": "SH-SY5Y\n(in vitro, human)",
}
PADJ_THRESH   = 0.05
FC_THRESH     = 0.5   # |log2FC| threshold for colouring

COL_UP     = "#d73027"
COL_DOWN   = "#4575b4"
COL_NS     = "#aaaaaa"
ALPHA_NS   = 0.08
ALPHA_SIG  = 0.55

# --- Type sizing -----------------------------------------------------------
# matplotlib sizes are absolute points on the canvas, but the manuscript scales
# the figure to the width it is placed at, so everything shrinks on the way to
# the page. That scale, not the numbers below, is what governs how large the
# text actually reads: on the original 13.5-inch canvas placed in one column it
# was 0.24, which rendered the 7 pt axis labels at 1.7 pt.
#
# So the sizes are declared as the points they should measure *in the printed
# manuscript* and the authored size is derived. Change PANEL_W and the type
# follows; change a target here and it means what it says.
#
# PLACED_WIDTH_IN must track how manuscript.tex includes the figure. It is a
# full-width float (figure* at \textwidth), because six volcano panels across
# one 8.2 cm column held every label under 5 pt however the canvas was sized.
PLACED_WIDTH_IN = 6.69      # 21 - 2*2 cm text width, from the geometry package
PANEL_W, PANEL_H = 3.2, 3.0
GRID_COLS = 3

#: Point sizes as they should appear on the manuscript page.
#
# The ceiling is the panel width, not taste. Three panels share the 482 pt text
# block, so each gets ~160 pt, and the longest title -- "(n=189: 96 case / 93
# control)", 29 characters -- needs roughly 3.4 pt per character before it runs
# into its neighbour. In a single column each panel had 77 pt and the same
# title capped every size here at 5 pt, which is why the figure went
# full-width.
ON_PAGE_PT = {
    "gene": 6.0,            # inline gene annotations, the smallest thing here
    "tick": 7.0,
    "axis": 7.5,
    "title": 8.0,
    "panel_letter": 10.0,
    "suptitle": 10.0,
}


def pt(element: str) -> float:
    """Authored point size that renders at ON_PAGE_PT[element] once placed."""
    return ON_PAGE_PT[element] * (GRID_COLS * PANEL_W) / PLACED_WIDTH_IN


def study_sample_size(study_id: str) -> str:
    """Return "n=A case / B control" for a study, or "n unavailable".

    The panels previously read n_case/n_control from effects_normalized_all.csv,
    which does not carry those columns, so every panel silently rendered
    "n=? samples". The per-study effect tables written by step 05 do carry
    them, so they are read from there instead.
    """
    path = PER_STUDY_DIR / f"{study_id}_effects.csv"
    if not path.exists():
        logger.warning("no per-study effects for %s; sample size unavailable", study_id)
        return "n unavailable"
    try:
        row = pd.read_csv(path, usecols=["n_case", "n_control"], nrows=1)
        n_case, n_control = int(row["n_case"][0]), int(row["n_control"][0])
    except (ValueError, KeyError, IndexError) as exc:
        logger.warning("could not read sample size for %s: %s", study_id, exc)
        return "n unavailable"
    return f"n={n_case + n_control}: {n_case} case / {n_control} control"


def volcano_ax(
    ax: plt.Axes,
    df: pd.DataFrame,
    title: str,
    fc_col: str = "effect_size",
    p_col: str  = "pval",
    padj_col: str = "padj",
    max_genes: int = 3,
) -> None:
    """Draw a single volcano panel."""
    df = df.dropna(subset=[fc_col, p_col]).copy()
    df["-log10p"] = -np.log10(df[p_col].clip(lower=1e-300))

    # Significance classification
    sig_up = (df[(df[padj_col] < PADJ_THRESH) & (df[fc_col] > FC_THRESH)]
              if padj_col in df.columns else pd.DataFrame())
    sig_down = (df[(df[padj_col] < PADJ_THRESH) & (df[fc_col] < -FC_THRESH)]
                if padj_col in df.columns else pd.DataFrame())
    ns       = df.drop(index=sig_up.index.tolist() + sig_down.index.tolist())

    # rasterized: these panels carry ~1.5M points in total, and one vector
    # path per point makes the figure unusable in any vector format (its SVG
    # rendition reached 469 MB). Only the point clouds are rasterized; axes,
    # ticks, labels and gene annotations stay vector, so all text remains
    # selectable and sharp at any zoom.
    ax.scatter(ns[fc_col], ns["-log10p"], s=2, alpha=ALPHA_NS, c=COL_NS,
               linewidths=0, zorder=1, rasterized=True)
    if not sig_up.empty:
        ax.scatter(sig_up[fc_col], sig_up["-log10p"], s=6, alpha=ALPHA_SIG,
                   c=COL_UP, linewidths=0, zorder=2, rasterized=True)
    if not sig_down.empty:
        ax.scatter(sig_down[fc_col], sig_down["-log10p"], s=6, alpha=ALPHA_SIG,
                   c=COL_DOWN, linewidths=0, zorder=2, rasterized=True)

    ax.axhline(-np.log10(0.05), color="grey", lw=0.6, ls="--", alpha=0.7)
    ax.axvline(FC_THRESH,  color=COL_UP,   lw=0.5, ls=":", alpha=0.5)
    ax.axvline(-FC_THRESH, color=COL_DOWN, lw=0.5, ls=":", alpha=0.5)

    # Label top genes by -log10p.
    #
    # The top hits of a volcano sit almost on top of one another, so labels
    # drawn at their own points overlap into an unreadable pile. They were too
    # small to notice before this figure's type was enlarged.
    #
    # Offsetting each label by its rank does not fix that: the offset is added
    # to a point that already has its own height, so a lower-ranked label can
    # still land on a higher-ranked one, which is what left three names bunched
    # in the pooled panel. The stack is therefore placed in axes coordinates --
    # fixed rows in a corner, one line apart by construction -- with a hairline
    # tying each name back to its own point. Nothing moves the points; only the
    # names are relocated, and each says which point it belongs to.
    #
    # The corners are free: set_ylim below leaves headroom above the tallest
    # point precisely so the stack has somewhere to sit.
    #
    # Row spacing is an offset in inches, not a fraction of the axes. Expressing
    # it as a fraction needs the axes height, and dividing by the panel height
    # instead -- which tight_layout shrinks the axes well below -- put the rows
    # about one line apart where 1.45 was intended, so they touched.
    line_in = pt("gene") * 1.45 / 72
    for sub, col in [(sig_up, COL_UP), (sig_down, COL_DOWN)]:
        if sub.empty:
            continue
        top = sub.nlargest(max_genes, "-log10p")
        # Up-regulated names in the top-right corner, down-regulated in the
        # top-left, matching the side of the plot their points are on.
        x, ha = (0.985, "right") if col == COL_UP else (0.015, "left")
        for rank, (_, row) in enumerate(top.iterrows()):
            stacked = offset_copy(ax.transAxes, fig=ax.figure,
                                  y=-rank * line_in, units="inches")
            ax.annotate(
                str(row["feature_id"]),
                xy=(row[fc_col], row["-log10p"]), xycoords="data",
                xytext=(x, 0.97), textcoords=stacked,
                fontsize=pt("gene"), color=col, ha=ha, va="top",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.3, alpha=0.55,
                                shrinkA=1, shrinkB=1),
                zorder=3,
            )

    n_up   = len(sig_up)
    n_down = len(sig_down)
    ax.set_title(f"{title}\n↑{n_up}  ↓{n_down}", fontsize=pt("title"), pad=4)
    ax.set_xlabel("log₂FC (case vs control)", fontsize=pt("axis"))
    ax.set_ylabel("-log₁₀(p)", fontsize=pt("axis"))
    ax.tick_params(labelsize=pt("tick"))

    # 1.18: margin for the gene labels, which start at their point and run
    # outward. Without it the outermost names are clipped mid-word by the axis.
    xmax = max(abs(df[fc_col]).quantile(0.999), FC_THRESH + 0.1) * 1.18
    ax.set_xlim(-xmax, xmax)
    # Room above the tallest point for the label ladder; without it the top
    # name is clipped by the axes.
    ymax = df["-log10p"].max()
    if np.isfinite(ymax) and ymax > 0:
        ax.set_ylim(top=ymax * (1.10 + 0.10 * max_genes))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    per_study = pd.read_csv(PER_STUDY_CSV)
    pooled    = pd.read_csv(POOLED_CSV)

    human_ps  = {sid: per_study[per_study["study_id"] == sid] for sid in HUMAN_STUDIES}

    n_panels = len(HUMAN_STUDIES) + 1   # +1 for rodent pooled

    cols = GRID_COLS
    rows = int(np.ceil(n_panels / cols))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(cols * PANEL_W, rows * PANEL_H))
    axes_flat = axes.flatten()

    panel_idx = 0

    # Pooled meta-analysis volcano. This is pooled_effects.csv unfiltered, so
    # it covers every contributing study, human ones included - it is not a
    # rodent-only panel. The study count is read from the data rather than
    # hardcoded, which is how the previous "n=21" survived the cohort
    # expansion to 59.
    ax = axes_flat[panel_idx]
    n_studies = per_study["study_id"].nunique()
    pooled_plot = pooled.rename(columns={"yi_pooled": "effect_size", "pval_pooled": "pval"})
    pooled_plot["padj"] = (pooled_plot["padj_pooled"]
                           if "padj_pooled" in pooled_plot.columns
                           else pooled_plot.get("padj", np.nan))
    volcano_ax(
        ax, pooled_plot,
        title=f"Pooled meta-analysis, all studies\n(k≥2, {n_studies} studies)",
        fc_col="effect_size", p_col="pval", padj_col="padj",
    )
    ax.set_facecolor("#f5f0e8")
    panel_idx += 1

    # Human study volcanoes
    for sid, label in HUMAN_STUDIES.items():
        if panel_idx >= len(axes_flat):
            break
        ax = axes_flat[panel_idx]
        df_h = human_ps.get(sid, pd.DataFrame())
        if df_h.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    fontsize=pt("axis"), transform=ax.transAxes)
            ax.set_title(label, fontsize=pt("title"))
        else:
            volcano_ax(
                ax, df_h,
                title=f"{label}\n({study_sample_size(sid)})",
                fc_col="effect_size", p_col="pval", padj_col="padj",
            )
        panel_idx += 1

    # Hide unused panels
    for ax in axes_flat[panel_idx:]:
        ax.set_visible(False)

    # Panel labels
    for i, ax in enumerate(axes_flat[:n_panels]):
        lbl = chr(ord("A") + i)
        ax.text(-0.17, 1.14, lbl, transform=ax.transAxes,
                fontsize=pt("panel_letter"), fontweight="bold", va="top")

    # Wrapped deliberately: as one line this measured wider than the canvas,
    # and bbox_inches="tight" then grew the figure to fit it, cutting the
    # placement scale from 0.336 to 0.189 and squeezing every panel. A title
    # must not be what decides the figure's width.
    fig.suptitle(
        "Transcriptomic differential expression:\n"
        "pooled meta-analysis and individual human studies",
        fontsize=pt("suptitle"), fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    # dpi governs only the rasterized point clouds here; 300 is the usual
    # journal floor for raster content, and 180 would be visibly soft in print.
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved → %s", OUT_PDF)

    # Print summary stats
    print("\nVolcano summary:")
    if "padj" in pooled_plot.columns:
        n_rod = (pooled_plot["padj"] < PADJ_THRESH).sum()
        print(f"  Rodent pooled: {n_rod} sig (padj<0.05)")
    for sid, label in HUMAN_STUDIES.items():
        df_h = human_ps.get(sid, pd.DataFrame())
        if df_h.empty or "padj" not in df_h.columns:
            continue
        n = (df_h["padj"] < PADJ_THRESH).sum()
        print(f"  {sid} ({label.split(chr(10))[0]}): {n} sig (padj<0.05), {len(df_h)} features")


if __name__ == "__main__":
    main()
