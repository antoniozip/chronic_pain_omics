"""PRISMA 2020 flow diagram, drawn per repository.

One row per source: identification on the left, exclusions with reasons in the
middle (branching off the main arrow, as PRISMA 2020 draws them), inclusion on
the right, and a closing box with what enters the synthesis. The counts come
from :mod:`cp_multiomics.prisma_flow`, which reads the screening sheets and
refuses a flow that does not add up.

The diagram is drawn at the journal's full print width with text at FONT, so
it is never scaled at print and every character prints at the size set here.
That is why each row shows only its largest exclusion reasons
(:meth:`SourceFlow.condensed`): the full breakdown could not fit a page at a
legible size, and it is released in the supplementary tables instead.

The figure is saved with a sidecar, ``prisma_flow_counts.csv``, holding every
number it draws. A test recomputes those numbers from the sheets, so a figure
left stale after a screening decision changes fails the suite instead of
reaching the journal: the figure itself is an image no claim can read.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath

from ..prisma_flow import MAX_REASONS, SourceFlow, flow_table
from .style import set_publication_style

logger = logging.getLogger(__name__)

# Every length below is in inches: the axes fill the figure and one data unit
# is one inch, so text and boxes keep the proportions they are laid out in.
# Under matplotlib's default margins the boxes shrank while the text did not.
#
# 11 pt is the largest size at which the diagram fits one page (about 229 mm
# of text height) at 174 mm wide. The journal recommends 12 pt, which made it
# 287 mm tall even with the reasons folded; 11 pt, 225 mm. The column widths
# are set by their widest words at FONT, "transcriptomics)" in bold on the
# left and "93 study units" on the right: change FONT and re-measure them,
# or _box raises on a word wider than its box.
FONT = 11.0
LINE = 1.13 * FONT / 72   # baseline to baseline
WIDTH = 174 / 25.4        # the journal's full print width
PAGE_HEIGHT = 229 / 25.4  # and about the text height of its page
MARGIN = 0.04             # keeps edges inside the saved page
PAD = 0.06                # between a box's edge and its text
INSET = 0.07
COUNT_GAP = 0.12          # between a label and its right-aligned count
COL_ID = (MARGIN, 1.66)
COL_EX = (1.84, 5.39)
COL_IN = (5.57, WIDTH - MARGIN)
ROW_GAP = 0.10
COLOURS = {"id": "#DCEAF7", "ex": "#F0F0F0", "await": "#FFF2CC", "in": "#C6DBEF"}
EDGE = "#2166AC"
_TEXT_TO_PATH = TextToPath()

# One logical line of a box: its label, the count set flush right beside it
# ("" for none), and whether it is set bold. Counts sit in their own column,
# as in a table, rather than as "(n = 98)" in the running text: that suffix
# alone pushed a dozen labels onto a second line.
Line = tuple[str, str, bool]


def _text_width(text: str, bold: bool = False) -> float:
    """Rendered width in inches, from the glyph outlines, on any backend."""
    prop = FontProperties(size=FONT, weight="bold" if bold else "normal")
    width, _, _ = _TEXT_TO_PATH.get_text_width_height_descent(text, prop, ismath=False)
    return width / 72


def _wrap(text: str, width_in: float, bold: bool = False) -> list[str]:
    """Break at spaces so every line measures at most `width_in` as rendered."""
    lines: list[str] = []
    for word in text.split(" "):
        if lines and _text_width(f"{lines[-1]} {word}", bold) <= width_in:
            lines[-1] += f" {word}"
        else:
            lines.append(word)
    return lines or [""]


def _box(ax: Axes, col: tuple[float, float], top: float, lines: list[Line],
         colour: str) -> float:
    """Draw a box hanging from `top`, as tall as its text; return its bottom."""
    x0, x1 = col
    counts = [_text_width(n, bold) for _, n, bold in lines if n]
    label_width = x1 - x0 - 2 * INSET - (max(counts) + COUNT_GAP if counts else 0.0)
    y = top - PAD
    texts = []
    for label, n, bold in lines:
        weight = "bold" if bold else "normal"
        if n:
            ax.text(x1 - INSET, y, n, ha="right", va="top", fontsize=FONT,
                    fontweight=weight)
        for part in _wrap(label, label_width, bold):
            # Only a single word can overrun, and it would print over the edge.
            if _text_width(part, bold) > label_width:
                raise ValueError(f"{part!r} is wider than its box")
            texts.append(ax.text(x0 + INSET, y, part, ha="left", va="top",
                                 fontsize=FONT, fontweight=weight))
            y -= LINE
    lowest = min(t.get_window_extent().y0 for t in texts) / ax.figure.dpi
    bottom = lowest - PAD
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, bottom), x1 - x0, top - bottom, boxstyle="round,pad=0,rounding_size=0.04",
        facecolor=colour, edgecolor=EDGE, linewidth=0.7, zorder=0))
    return bottom


def _items(items: tuple[tuple[str, int], ...]) -> list[Line]:
    return [(label, f"{n:,}", False) for label, n in items]


def _arrow(ax: Axes, x0: float, y0: float, x1: float, y1: float) -> None:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.8,
                                shrinkA=0, shrinkB=0))


def _row(ax: Axes, flow: SourceFlow, top: float) -> float:
    """Draw one source's row; return the lowest y it used."""
    id_lines: list[Line] = [(flow.source, "", True),
                            (f"{flow.identified:,} {flow.identified_label}", "", False)]
    in_lines: list[Line] = [(f"Included: {flow.included:,} {flow.included_label}", "", False)]
    if flow.units is not None:
        in_lines.append((f"{flow.units:,} study units", "", False))
    bottoms = [_box(ax, COL_ID, top, id_lines, COLOURS["id"]),
               _box(ax, COL_IN, top, in_lines, COLOURS["in"])]

    # The main arrow runs at the height of the first text line; the exclusions
    # branch down from it into the middle column.
    y_arrow = top - PAD - 0.45 * FONT / 72
    _arrow(ax, COL_ID[1], y_arrow, COL_IN[0], y_arrow)
    x_branch = 0.5 * (COL_EX[0] + COL_EX[1])
    ex_top = y_arrow - 0.08
    _arrow(ax, x_branch, y_arrow, x_branch, ex_top)

    items = flow.removed + flow.excluded
    title = "Removed or excluded" if flow.removed else "Excluded"
    y = _box(ax, COL_EX, ex_top,
             [(title, f"{sum(n for _, n in items):,}", True)] + _items(items),
             COLOURS["ex"])
    if flow.awaiting:
        y = _box(ax, COL_EX, y - 0.06,
                 [("Awaiting classification", f"{flow.n_awaiting:,}", True)]
                 + _items(flow.awaiting), COLOURS["await"])
    bottoms.append(y)
    return min(bottoms)


def plot_repository_flow(flows: list[SourceFlow], summary: list[str],
                         out_dir: Path, max_reasons: int = MAX_REASONS) -> Path:
    """Draw the per-repository flow diagram and its counts sidecar.

    Args:
        flows: One :class:`SourceFlow` per source, in drawing order.
        summary: Lines for the closing "included in the synthesis" box.
        out_dir: Output directory.
        max_reasons: Exclusion reasons drawn per source before the rest fold
            into one line; see :meth:`SourceFlow.condensed`.

    Returns:
        Path to the saved PDF.
    """
    set_publication_style()
    drawn = [flow.condensed(max_reasons) for flow in flows]
    for flow in drawn:
        flow.check()

    # Lay out top-down from a generous height, then crop the figure to what
    # was drawn, keeping one data unit to one inch.
    height = 20.0
    fig = plt.figure(figsize=(WIDTH, height))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, height)
    ax.axis("off")

    top = height - MARGIN
    for col, label in ((COL_ID, "Identification"),
                       (COL_EX, "Screening and eligibility"),
                       (COL_IN, "Included")):
        ax.text(0.5 * sum(col), top, label, ha="center", va="top", fontsize=FONT,
                fontweight="bold", color=EDGE)
    top -= 1.45 * FONT / 72
    for flow in drawn:
        top = _row(ax, flow, top) - ROW_GAP

    lines: list[Line] = [("Included in the synthesis", "", True)]
    lines += [(item, "", False) for item in summary]
    bottom = _box(ax, (MARGIN, WIDTH - MARGIN), top, lines, COLOURS["in"]) - MARGIN

    ax.set_ylim(bottom, height)
    fig.set_size_inches(WIDTH, height - bottom)
    if height - bottom > PAGE_HEIGHT:
        logger.warning("PRISMA flow is %.0f mm tall, over a page: it will be scaled "
                       "at print and its text will print below %g pt",
                       (height - bottom) * 25.4, FONT)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "prisma_flow.pdf"
    # Saved at the figure's own size, not the style's tight crop, which pads
    # the drawing and would make the page wider than the print width.
    fig.savefig(out_path, bbox_inches=fig.bbox_inches)
    plt.close(fig)
    flow_table(drawn).to_csv(out_dir / "prisma_flow_counts.csv", index=False)
    logger.info("Saved PRISMA flow: %s (%.0f x %.0f mm)", out_path,
                WIDTH * 25.4, (height - bottom) * 25.4)
    return out_path
