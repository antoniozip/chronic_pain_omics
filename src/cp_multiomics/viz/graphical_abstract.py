"""Graphical abstract: the three disjointness axes named in the title.

The title claims that pain models, species and genomic risk implicate
non-overlapping gene sets, so the figure shows those three axes and the corpus
they rest on, rather than illustrating the pipeline. Every number is read from
``results/`` at draw time and none is passed in, so the figure cannot drift
from the text the way a hand-drawn asset would. A panel whose source file is
missing stops the whole figure rather than leaving a gap that would read as a
finding.

Panel A is the corpus. The genomic arm is counted differently from the others
by necessity: it has no pooled effect table because it is a replication-count
meta-analysis over gene-level statistics, so its bar is a gene count while the
others are correction families. The axis label says so, because a bar that
means something different from its neighbours and does not admit it is a
figure that misleads.

Panel B is model specificity. SNL dominates, and the Results explain why: it is
the best-covered stratum, so this is depth rather than biology.

Panel C is cross-species concordance against the 50% line. Both species sit
just below chance.

The genomic-transcriptomic disjointness, the title's third axis, is *not* drawn.
`results/meta/genomics/gwas_tx_overlap.csv` counts GWAS genes present in each
measured transcript set, which is coverage rather than agreement, and reads
59 of 142 for human transcripts -- the opposite of the paper's claim if shown
without that distinction. The claim itself rests on the 25 MAGMA-significant
genes, whose Entrez identifiers cannot be joined to the symbol-keyed
transcriptomic pool without a mapping this repository does not carry. It is
stated in the Results with proper support and left out here rather than
approximated.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

from ..strata import plain_label
from .style import set_publication_style

logger = logging.getLogger(__name__)

# Modality accent colours, kept distinct in greyscale as well as in colour:
# a printed reader and a screen reader should both be able to tell the arms
# apart.
ARM_COLOURS = {
    "Transcriptomic": "#2c5f8a",
    "Genomic": "#7a9cc6",
    "Proteomic": "#b8654f",
    "Metabolomic": "#5b8c5a",
}
STRATUM_COLOUR = "#2c5f8a"
CHANCE_COLOUR = "#b8654f"
MUTED = "#555555"


def _family_and_significant(path: Path) -> tuple[int, int]:
    """(features in the correction family, significant at padj < 0.05).

    The family is the set carrying an adjusted p-value at all. Features below
    the study threshold are pass-throughs of a single study's raw p-value and
    are deliberately excluded from correction, so counting them here would
    inflate every arm.
    """
    df = pd.read_csv(path, usecols=lambda c: c in {"padj_pooled", "padj"})
    col = "padj_pooled" if "padj_pooled" in df.columns else "padj"
    return int(df[col].notna().sum()), int((df[col] < 0.05).sum())


def _genomic_gene_count(meta_dir: Path) -> tuple[int, int]:
    """(genes carried forward, genes replicated in >1 study).

    The genomic arm never produced pooled effect sizes: no study in the
    retrieved set released usable summary statistics for a common variant set,
    so it is a replication-count meta-analysis over gene-level results. Its bar
    therefore counts genes, not a correction family.
    """
    df = pd.read_csv(meta_dir / "genomics" / "gwas_gene_summary.csv",
                     usecols=["gene", "n_studies"])
    return len(df), int((df["n_studies"] > 1).sum())


def _panel_corpus(ax, meta_dir: Path, arms: list[tuple[str, str, int]]) -> None:
    labels, families, sigs, colours = [], [], [], []
    for label, modality, n_studies in arms:
        if modality == "genomics":
            family, sig = _genomic_gene_count(meta_dir)
        else:
            path = meta_dir / modality / "pooled_effects.csv"
            if not path.exists():
                raise FileNotFoundError(path)
            family, sig = _family_and_significant(path)
        labels.append(f"{label}\n{n_studies} studies")
        families.append(max(family, 1))       # log axis cannot show 0
        sigs.append(sig)
        colours.append(ARM_COLOURS[label])

    y = range(len(labels))
    ax.barh(list(y), families, color=colours, height=0.62, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("features in the correction family, or genes for the\n"
                  "genomic arm (log scale)")
    ax.set_title("A  Evidence assembled", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.9", zorder=0)
    ax.set_axisbelow(True)

    # The genomic arm's second number is genes replicated in more than one
    # study, not significant features: it has no pooled effect sizes, so
    # writing "sig." beside it would claim a test that was never run.
    for i, (label, fam, sig) in enumerate(zip(labels, families, sigs, strict=True)):
        note = f"{sig:,} replicated" if label.startswith("Genomic") else f"{sig:,} sig."
        ax.text(fam * 1.35, i, f"{fam:,}  ({note})", va="center",
                fontsize=7.5, color=MUTED)
    ax.set_xlim(right=max(families) * 12)


def _panel_models(ax, strat_dir: Path) -> None:
    rows = []
    for path in sorted(strat_dir.glob("*_pooled.csv")):
        name = path.stem.replace("_pooled", "")
        if name.endswith("_human") or name == "invitro_human":
            continue                      # single-study strata; no family
        df = pd.read_csv(path, usecols=["padj"])
        rows.append((name, int((df["padj"] < 0.05).sum())))
    rows.sort(key=lambda r: -r[1])

    # The same words the stratified table prints, so a reader comparing the
    # two is not left matching "other_rodent_amendment" against
    # "Other (rodent)" by eye.
    names = [plain_label(r[0]) for r in rows]
    counts = [r[1] for r in rows]
    bars = ax.bar(names, [max(c, 0.4) for c in counts],
                  color=STRATUM_COLOUR, width=0.62, zorder=3)
    ax.set_yscale("log")
    ax.set_ylabel("significant features")
    ax.set_title("B  Pain models are near-disjoint", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0.3, max(counts) * 6)
    # Stratum names are long enough to overlap at this panel's width; rotating
    # them upright keeps each on one line instead of truncating or wrapping.
    ax.tick_params(axis="x", labelrotation=90)

    for bar, count in zip(bars, counts, strict=True):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.25,
                f"{count:,}", ha="center", fontsize=7.5, color=MUTED)


def _panel_species(ax, summary: Path) -> None:
    df = pd.read_csv(summary)
    df = df.sort_values("species")
    names = [{"Mus musculus": "Mouse", "Rattus norvegicus": "Rat"}.get(s, s)
             for s in df["species"]]
    pct = df["pct_concordant"].astype(float).tolist()

    bars = ax.bar(names, pct, color=STRATUM_COLOUR, width=0.5, zorder=3)
    ax.axhline(50, color=CHANCE_COLOUR, linestyle="--", linewidth=1.1, zorder=4)
    ax.text(ax.get_xlim()[1], 50, "  chance", va="center", ha="left",
            fontsize=7.5, color=CHANCE_COLOUR)
    ax.set_ylim(40, 53)
    ax.set_ylabel("directional concordance (%)")
    ax.set_title("C  Species agree at chance", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", zorder=0)
    ax.set_axisbelow(True)

    for bar, value, n in zip(bars, pct, df["n_ortholog_pairs"], strict=True):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.45,
                f"{value:.1f}%", ha="center", fontsize=8, fontweight="bold")
        # Inside the bar: at this panel's width the labels collide if set
        # below the axis, and an ortholog count is context rather than a
        # reading, so it should not compete with the percentage.
        ax.text(bar.get_x() + bar.get_width() / 2, 41.2,
                f"{int(n):,}\northologs", ha="center", va="bottom",
                fontsize=6.8, color="white", zorder=5, linespacing=1.25)


def plot_graphical_abstract(results_dir: Path, out_path: Path,
                            arms: list[tuple[str, str, int]]) -> Path | None:
    """Draw the graphical abstract, or return None if an input is missing."""
    set_publication_style()
    meta = results_dir / "meta"
    strat = meta / "transcriptomics" / "stratified"
    species = results_dir / "cross_species" / "transcriptomics" / "summary.csv"
    for required in (meta, strat, species):
        if not required.exists():
            logger.warning("graphical abstract: %s missing; not drawn", required)
            return None

    # Panel B's stratum names are rotated upright, which costs roughly an inch
    # of height. Without it the panels are squeezed into the top third and
    # panel C's y-label runs past the axis it belongs to.
    fig_h = 4.35
    fig, axes = plt.subplots(1, 3, figsize=(10.6, fig_h),
                             gridspec_kw={"width_ratios": [1.32, 1.0, 0.80]})
    _panel_corpus(axes[0], meta, arms)
    _panel_models(axes[1], strat)
    _panel_species(axes[2], species)

    # A single rule under the panels ties them into one statement rather than
    # three unrelated charts. The caption band is a fixed height in inches, so
    # it does not grow with the figure.
    band = 0.39 / fig_h
    fig.add_artist(Rectangle((0.045, band * 0.63), 0.915, 0.004,
                             transform=fig.transFigure, color="0.82", lw=0))
    fig.text(0.045, band * 0.13,
             "Evidence layers in chronic pain do not converge: pain models share "
             "few significant features, and rodent and human\ntranscriptomes agree "
             "no better than chance.",
             fontsize=8.5, color="0.25", transform=fig.transFigure)

    fig.tight_layout(rect=(0, band, 1, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    logger.info("graphical abstract -> %s", out_path)
    return out_path
