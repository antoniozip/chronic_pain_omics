"""Publication-style matplotlib rcParams for all figures."""

from __future__ import annotations

import matplotlib as mpl
from cycler import cycler


def set_publication_style() -> None:
    """Apply consistent publication-quality style to all subsequent plots."""
    mpl.rcParams.update(
        {
            # Font
            "font.family":       "sans-serif",
            "font.sans-serif":   ["Arial", "DejaVu Sans"],
            "font.size":         9,
            "axes.titlesize":    10,
            "axes.labelsize":    9,
            "xtick.labelsize":   8,
            "ytick.labelsize":   8,
            "legend.fontsize":   8,
            # Figure
            "figure.dpi":        300,
            "savefig.dpi":       300,
            "savefig.format":    "pdf",
            "savefig.bbox":      "tight",
            # Lines / spines
            "axes.linewidth":    0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "lines.linewidth":   1.2,
            # Grid
            "axes.grid":         False,
            # Colors
            "axes.prop_cycle":   cycler(
                color=["#2166AC", "#D6604D", "#4DAC26", "#8073AC", "#E08214"]
            ),
        }
    )


# Significance color palette (shared by the volcano and concordance plots)
SIG_COLORS = {
    "ns":               "#AAAAAA",
    "sig":              "#2166AC",
    "sig_up":           "#D6604D",
    "sig_down":         "#2166AC",
    "both_sig_up":      "#B2182B",
    "both_sig_down":    "#053061",
}
