"""GWAS-transcriptomics figure: replicated risk genes and what the pool says.

Generates:
  manuscript/figures/gwas_tx_overlap.pdf
    (A) the most-replicated GWAS pain genes, shaded by independent studies
    (B) the pooled transcriptomic effect of those same genes

Replaces the draft figure of the same name, whose generator was never
committed and whose content had gone stale: it was built on the May 2026
transcriptomic pool, in which these genes carried k=3, and the August pool
gives every one of them a different effect size at k=15-23, ASTN2 with the
opposite sign. Its gene panel also plotted ten genes in alphabetical order
under the label "overlap genes" by a rule nothing in the repository records.

Both panels here are drawn from the tracked results, and both are checked by
tests/test_manuscript_claims.py so this cannot happen quietly again.

The two panels share one gene set on purpose. The paper's claim is that the
genomic and transcriptomic layers implicate disjoint genes, and the honest way
to show it is to take the genes genetics nominates most strongly and ask what
the transcriptomic meta-analysis says about exactly those -- not to plot a
subset chosen after the fact.

    uv run python pipeline/15_gwas_tx_overlap.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
REPLICATION_CSV = REPO_ROOT / "results" / "meta" / "genomics" / "gwas_gene_meta_replication.csv"
POOLED_CSV = REPO_ROOT / "results" / "meta" / "transcriptomics" / "pooled_effects.csv"
OUT_PDF = REPO_ROOT / "manuscript" / "figures" / "gwas_tx_overlap.pdf"

#: Genes shown. The replication table runs to 142; the figure shows the head of
#: it, which is where the genetic evidence is strong enough for the absence of
#: a transcriptomic counterpart to mean anything.
TOP_N = 15
PADJ_THRESH = 0.05

COL_UP = "#d73027"
COL_DOWN = "#4575b4"

# Sizes are the points they should measure on the manuscript page; the authored
# size follows from the width the figure is placed at, as in
# pipeline/13_human_rodent_volcano.py. This figure is a single-column float.
PLACED_WIDTH_IN = 3.23
FIG_W, FIG_H = 6.0, 6.6
ON_PAGE_PT = {"tick": 5.5, "axis": 6.0, "title": 6.5, "panel_letter": 8.5}


def pt(element: str) -> float:
    """Authored point size rendering at ON_PAGE_PT[element] once placed."""
    return ON_PAGE_PT[element] * FIG_W / PLACED_WIDTH_IN


def load_genes(top_n: int = TOP_N) -> pd.DataFrame:
    """The most-replicated GWAS pain genes, ordered by replication count.

    Ties are broken on the minimum p-value, so the order is defined rather
    than dependent on the row order the table happens to carry.
    """
    rep = pd.read_csv(REPLICATION_CSV, usecols=["gene", "n_studies", "min_p"])
    rep = rep.sort_values(["n_studies", "min_p"], ascending=[False, True])
    return rep.head(top_n).reset_index(drop=True)


def load_pooled_effects(genes: list[str]) -> pd.DataFrame:
    """Pooled transcriptomic estimate for each gene, NaN where not pooled."""
    pooled = pd.read_csv(POOLED_CSV,
                         usecols=["feature_id", "k", "yi_pooled", "padj_pooled"])
    pooled = pooled[pooled["feature_id"].isin(genes)]
    return pooled.set_index("feature_id").reindex(genes)


def _panel_replication(ax: plt.Axes, rep: pd.DataFrame) -> None:
    order = rep.iloc[::-1]                       # largest at the top
    norm = Normalize(vmin=order["n_studies"].min(), vmax=order["n_studies"].max())
    colours = plt.cm.Blues(norm(order["n_studies"]))
    ax.barh(order["gene"], order["n_studies"], color=colours,
            edgecolor="#333333", linewidth=0.4)
    ax.set_xlabel("Independent GWAS reporting the gene", fontsize=pt("axis"))
    ax.set_title("Most replicated GWAS pain genes", fontsize=pt("title"), pad=4)
    ax.tick_params(labelsize=pt("tick"))
    ax.set_xlim(0, order["n_studies"].max() * 1.08)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _panel_transcriptomic(ax: plt.Axes, rep: pd.DataFrame,
                          eff: pd.DataFrame) -> int:
    """Pooled log2FC of the same genes. Returns how many are missing."""
    order = rep.iloc[::-1]
    values = eff.loc[order["gene"], "yi_pooled"]
    absent = int(values.isna().sum())
    plotted = values.fillna(0.0)
    colours = [COL_UP if v > 0 else COL_DOWN for v in plotted]
    ax.barh(order["gene"], plotted, color=colours, edgecolor="none")
    ax.axvline(0, color="#333333", lw=0.6)
    ax.set_xlabel("log$_2$FC (pooled transcriptomic meta-analysis)",
                  fontsize=pt("axis"))
    ax.set_title("Transcriptomic effect of the same genes",
                 fontsize=pt("title"), pad=4)
    ax.tick_params(labelsize=pt("tick"))
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    return absent


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    rep = load_genes()
    eff = load_pooled_effects(list(rep["gene"]))

    fig, axes = plt.subplots(2, 1, figsize=(FIG_W, FIG_H))
    _panel_replication(axes[0], rep)
    absent = _panel_transcriptomic(axes[1], rep, eff)

    for letter, ax in zip("AB", axes):
        ax.text(-0.28, 1.08, letter, transform=ax.transAxes,
                fontsize=pt("panel_letter"), fontweight="bold", va="top")

    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved → %s", OUT_PDF)

    sig = eff["padj_pooled"] < PADJ_THRESH
    logger.info("%d genes shown; %d absent from the pooled table; "
                "%d significant at padj<%.2f",
                len(rep), absent, int(sig.sum()), PADJ_THRESH)
    logger.info("pooled k: %s",
                ", ".join(f"{g}={'-' if np.isnan(k) else int(k)}"
                          for g, k in eff["k"].items()))


if __name__ == "__main__":
    main()
