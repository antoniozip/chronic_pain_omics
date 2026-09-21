"""Visualization package for paper-quality figures.

All functions write to manuscript/figures/ and return the output Path.
Global rcParams are set once via set_publication_style().
"""

from .concordance_plot import plot_concordance_scatter
from .graphical_abstract import plot_graphical_abstract
from .heatmap import plot_cross_modal_heatmap
from .prisma import plot_repository_flow
from .style import set_publication_style
from .volcano import plot_volcano

__all__ = [
    "set_publication_style",
    "plot_volcano",
    "plot_cross_modal_heatmap",
    "plot_concordance_scatter",
    "plot_repository_flow",
    "plot_graphical_abstract",
]
