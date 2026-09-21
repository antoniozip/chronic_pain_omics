"""PRISMA 2020 flow diagram, drawn per repository.

One row per source: identification on the left, exclusions with reasons in the
middle (branching off the main arrow, as PRISMA 2020 draws them), inclusion on
the right, and a closing box with what enters the synthesis. The counts come
from :mod:`cp_multiomics.prisma_flow`, which reads the screening sheets and
refuses a flow that does not add up.

The figure is saved with a sidecar, ``prisma_flow_counts.csv``, holding every
number it draws. A test recomputes those numbers from the sheets, so a figure
left stale after a screening decision changes fails the suite instead of
reaching the journal: the figure itself is an image no claim can read.
"""

from __future__ import annotations

import logging
import textwrap
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from ..prisma_flow import SourceFlow, flow_table
from .style import set_publication_style

logger = logging.getLogger(__name__)

FONT = 8.0
LINE = 0.148          # inches per text line at FONT
PAD = 0.12
WIDTH = 7.5
COL_ID = (0.0, 2.0)
COL_EX = (2.25, 5.2)
COL_IN = (5.6, WIDTH)
COLOURS = {"id": "#DCEAF7", "ex": "#F0F0F0", "await": "#FFF2CC", "in": "#C6DBEF"}
EDGE = "#2166AC"


def _wrap(text: str, width_in: float) -> list[str]:
    chars = max(12, int(width_in * 72 / (FONT * 0.60)))
    return textwrap.wrap(text, chars) or [""]


def _lines(title: str, items: tuple[tuple[str, int], ...], width_in: float) -> list[str]:
    out = [title]
    for label, n in items:
        # Non-breaking spaces keep "(n = 98)" whole; wrapped alone it reads as
        # a separate item.
        out += _wrap(f"{label} (n\u00a0=\u00a0{n:,})", width_in - 0.2)
    return out


def _box(ax, x0: float, x1: float, top: float, lines: list[str], colour: str) -> float:
    """Draw a box hanging from `top`; return its bottom."""
    h = len(lines) * LINE + 2 * PAD
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, top - h), x1 - x0, h, boxstyle="round,pad=0.02",
        facecolor=colour, edgecolor=EDGE, linewidth=0.7))
    ax.text(x0 + 0.08, top - PAD, "\n".join(lines), ha="left", va="top",
            fontsize=FONT, linespacing=1.25)
    return top - h


def _arrow(ax, x0: float, y0: float, x1: float, y1: float) -> None:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.8))


def _row(ax, flow: SourceFlow, top: float) -> float:
    """Draw one source's row; return the lowest y it used."""
    id_lines = _wrap(flow.source, COL_ID[1] - COL_ID[0] - 0.15) + [
        f"{flow.identified:,} {flow.identified_label}"]
    in_lines = [f"Included (n = {flow.included:,})", flow.included_label]
    if flow.units is not None:
        in_lines.append(f"{flow.units:,} study units")
    bottoms = [_box(ax, *COL_ID, top, id_lines, COLOURS["id"]),
               _box(ax, *COL_IN, top, in_lines, COLOURS["in"])]

    y_arrow = top - PAD - 0.5 * LINE
    _arrow(ax, COL_ID[1] + 0.05, y_arrow, COL_IN[0] - 0.05, y_arrow)
    x_branch = 0.5 * (COL_EX[0] + COL_EX[1])
    ex_top = y_arrow - 0.18
    _arrow(ax, x_branch, y_arrow, x_branch, ex_top + 0.02)

    width = COL_EX[1] - COL_EX[0]
    items = flow.removed + flow.excluded
    title = f"Removed or excluded (n = {sum(n for _, n in items):,})"
    y = _box(ax, *COL_EX, ex_top, _lines(title, items, width), COLOURS["ex"])
    if flow.awaiting:
        y = _box(ax, *COL_EX, y - 0.08,
                 _lines(f"Awaiting classification (n = {flow.n_awaiting:,})",
                        flow.awaiting, width), COLOURS["await"])
    bottoms.append(y)
    return min(bottoms)


def plot_repository_flow(flows: list[SourceFlow], summary: list[str],
                         out_dir: Path) -> Path:
    """Draw the per-repository flow diagram and its counts sidecar.

    Args:
        flows: One :class:`SourceFlow` per source, in drawing order.
        summary: Lines for the closing "included in the synthesis" box.
        out_dir: Output directory.

    Returns:
        Path to the saved PDF.
    """
    set_publication_style()
    # Lay out top-down in inch units, then size the figure to what was drawn.
    fig, ax = plt.subplots(figsize=(WIDTH, 11.0))
    height = 11.0
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, height)
    ax.axis("off")

    top = height - 0.05
    for x, label in ((0.5 * sum(COL_ID), "Identification"),
                     (0.5 * sum(COL_EX), "Screening and eligibility"),
                     (0.5 * sum(COL_IN), "Included")):
        ax.text(x, top, label, ha="center", va="top", fontsize=FONT + 1,
                fontweight="bold", color=EDGE)
    top -= 0.35
    for flow in flows:
        top = _row(ax, flow, top) - 0.22

    lines: list[str] = ["Included in the synthesis"]
    for item in summary:
        lines += _wrap(item, WIDTH - 0.3)
    bottom = _box(ax, 0.0, WIDTH, top, lines, COLOURS["in"])

    ax.set_ylim(bottom - 0.05, height)
    fig.set_size_inches(WIDTH, height - bottom + 0.05)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "prisma_flow.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    flow_table(flows).to_csv(out_dir / "prisma_flow_counts.csv", index=False)
    logger.info("Saved PRISMA flow: %s", out_path)
    return out_path
