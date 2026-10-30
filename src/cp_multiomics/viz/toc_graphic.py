"""Table-of-contents image: rodent-human agreement against a within-species floor.

iMetaOmics prints a table-of-contents entry beside each abstract, with an
image that must "fit within the dimensions of 50mm x 60mm and be fully legible
at this size". At that size one comparison is all that reads, so the image
shows the paper's central control rather than summarising every arm: how often
two pooled analyses agree on the direction of change, for two halves of one
human condition and for each rodent species against the human pool, with chance
marked.

Every value is read from ``results/`` at draw time, and the values plotted are
written beside the image (``toc_graphic_values.csv``) because a picture is not
something the claim registry can read; ``test_toc_graphic_plots_the_claimed_values``
binds that file to the claims stating the same quantities in the abstract.

The canvas is drawn at exactly 50 x 60 mm. The publication style saves with a
tight bounding box, which would crop the canvas to its contents and lose the
size the journal asks for, so this figure saves with the standard box.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from .style import set_publication_style

logger = logging.getLogger(__name__)

MM = 1 / 25.4
WIDTH_MM, HEIGHT_MM = 50, 60
CHANCE = 50.0

# The colours of the graphical abstract: blue for the comparison, the muted
# red-brown for chance. No red against green.
POINT_COLOUR = "#2c5f8a"
CHANCE_COLOUR = "#b8654f"
MUTED = "#555555"

SPECIES = {"mouse": "Mus musculus", "rat": "Rattus norvegicus"}


def _pct_concordant(summary: Path, species: str) -> float:
    d = pd.read_csv(summary)
    row = d[d["species"] == species]
    if row.empty:
        raise ValueError(f"{summary} has no row for {species}")
    return float(row["pct_concordant"].iloc[0])


def toc_values(cross_species_dir: Path) -> pd.Series:
    """The five percentages the image plots, read from the results.

    ``human_halves`` is the mean over the splits of the endometriosis stratum,
    the design the abstract quotes; the rodent values come from the whole
    human pool and from the tissue-matched sensitivity arm.
    """
    hvh = pd.read_csv(cross_species_dir / "human_vs_human.csv")
    replicate = hvh.loc[hvh["design"] == "replicate", "pct_concordant"]
    if replicate.empty:
        raise ValueError("human_vs_human.csv has no replicate design")
    all_human = cross_species_dir / "summary.csv"
    matched = cross_species_dir / "sensitivity" / "tissue_matched" / "summary.csv"
    values = {"human_halves": float(replicate.mean())}
    for name, species in SPECIES.items():
        values[f"{name}_all_human"] = _pct_concordant(all_human, species)
        values[f"{name}_tissue_matched"] = _pct_concordant(matched, species)
    return pd.Series(values, name="value")


def plot_toc_graphic(results_dir: Path, out_path: Path) -> Path | None:
    """Draw the image and its values sidecar, or return None if an input is missing."""
    cs_dir = results_dir / "cross_species" / "transcriptomics"
    for required in (cs_dir / "human_vs_human.csv", cs_dir / "summary.csv",
                     cs_dir / "sensitivity" / "tissue_matched" / "summary.csv"):
        if not required.exists():
            logger.warning("table-of-contents image: %s missing; not drawn", required)
            return None
    v = toc_values(cs_dir)

    set_publication_style()
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM))
    ax = fig.add_axes((0.27, 0.30, 0.67, 0.50))

    # Two-line labels: on one line "Mouse vs human" is wider than the margin
    # a 50 mm canvas can give it.
    rows = [("Human\nhalves", ["human_halves"]),
            ("Rat vs\nhuman", ["rat_all_human", "rat_tissue_matched"]),
            ("Mouse vs\nhuman", ["mouse_all_human", "mouse_tissue_matched"])]
    for y, (_, series) in enumerate(reversed(rows)):
        for s in series:
            filled = not s.endswith("tissue_matched")
            ax.plot(v[s], y, "o", ms=4.2, mew=0.9, color=POINT_COLOUR,
                    mfc=POINT_COLOUR if filled else "white", zorder=3)
    ax.set_yticks(range(len(rows)), [label for label, _ in reversed(rows)],
                  fontsize=6.5, linespacing=1.1)
    ax.set_ylim(-0.6, len(rows) - 0.4)

    ax.axvline(CHANCE, color=CHANCE_COLOUR, ls="--", lw=0.8, zorder=1)
    ax.text(CHANCE, len(rows) - 0.35, "chance", color=CHANCE_COLOUR,
            fontsize=6, ha="center", va="bottom")
    ax.axvline(v["human_halves"], color=MUTED, ls=":", lw=0.8, zorder=1)
    ax.text(v["human_halves"], len(rows) - 0.35, "floor", color=MUTED,
            fontsize=6, ha="center", va="bottom")

    lo = min(v.min(), CHANCE) - 1.5
    hi = max(v.max(), CHANCE) + 1.5
    ax.set_xlim(lo, hi)
    ax.tick_params(axis="x", labelsize=6, length=2, width=0.6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Same direction of change (%)", fontsize=6.5, labelpad=2)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)

    fig.text(0.5, 0.965, "Rodent\u2013human agreement\nfalls below a human floor",
             ha="center", va="top", fontsize=7, fontweight="bold",
             linespacing=1.15)
    # The key sits under the axis rather than in a legend box, which at this
    # size would cover a point.
    fig.text(0.06, 0.075, "● all human studies\n○ tissue-matched human studies",
             ha="left", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.3)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(out_path)
        fig.savefig(out_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
    v.rename_axis("series").reset_index().to_csv(
        out_path.with_name(out_path.stem + "_values.csv"), index=False)
    logger.info("table-of-contents image -> %s", out_path)
    return out_path
