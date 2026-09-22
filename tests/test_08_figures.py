"""Tests for viz package and pipeline/08_figures.py."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")  # non-interactive backend for tests

from cp_multiomics.viz.concordance_plot import plot_concordance_scatter
from cp_multiomics.viz.style import SIG_COLORS, set_publication_style
from cp_multiomics.viz.volcano import plot_volcano

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _pooled_df(n=50, seed=42):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "feature_id":   [f"GENE{i}" for i in range(n)],
        "yi_pooled":    rng.normal(0, 1, n),
        "padj_pooled":  np.concatenate([rng.uniform(0, 0.04, 10),
                                         rng.uniform(0.05, 1.0, n - 10)]),
        "k":            rng.integers(2, 10, n),
        "I2":           rng.uniform(0, 100, n),
    })


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def test_set_publication_style_runs():
    set_publication_style()  # should not raise
    import matplotlib as mpl
    assert mpl.rcParams["figure.dpi"] == 300


def test_sig_colors_defined():
    assert "sig_up" in SIG_COLORS
    assert "sig_down" in SIG_COLORS
    assert "ns" in SIG_COLORS


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

def test_plot_volcano_creates_pdf(tmp_path):
    df = _pooled_df()
    out = plot_volcano(df, modality="transcriptomics", out_dir=tmp_path,
                       padj_threshold=0.05, lfc_threshold=0.5)
    assert out.exists()
    assert out.suffix == ".pdf"


def test_plot_volcano_handles_empty_df(tmp_path):
    df = pd.DataFrame(columns=["feature_id", "yi_pooled", "padj_pooled"])
    out = plot_volcano(df, modality="transcriptomics", out_dir=tmp_path)
    # Should still produce a file (empty scatter)
    assert out.exists()


def test_plot_volcano_all_ns(tmp_path):
    df = _pooled_df()
    df["padj_pooled"] = 1.0   # nothing significant
    out = plot_volcano(df, modality="proteomics", out_dir=tmp_path)
    assert out.exists()


# ---------------------------------------------------------------------------
# Concordance scatter
# ---------------------------------------------------------------------------

def test_plot_concordance_scatter_creates_pdf(tmp_path):
    df = pd.DataFrame({
        "human_feature_id": [f"GENE{i}" for i in range(20)],
        "human_yi":         np.random.normal(0, 1, 20),
        "animal_yi":        np.random.normal(0, 1, 20),
        "human_padj":       np.random.uniform(0, 1, 20),
        "animal_padj":      np.random.uniform(0, 1, 20),
        "concordant":       np.random.choice([True, False], 20),
        "both_significant": np.random.choice([True, False], 20),
    })
    out = plot_concordance_scatter(
        df, modality="transcriptomics",
        animal_species="Mus musculus", out_dir=tmp_path,
    )
    assert out.exists()


def test_plot_concordance_scatter_empty_df(tmp_path):
    out = plot_concordance_scatter(
        pd.DataFrame(), modality="transcriptomics",
        animal_species="Mus musculus", out_dir=tmp_path,
    )
    # Returns a sentinel path, doesn't raise
    assert out is not None


# ---------------------------------------------------------------------------
# PRISMA flow, per repository
# ---------------------------------------------------------------------------

REPO = Path(__file__).resolve().parent.parent


def test_repository_flows_add_up_and_match_the_drawn_figure():
    """Every source adds up, and the figure states what the sheets say now.

    The figure is an image, so no claim can read it. Its sidecar records every
    number drawn; recomputing them here fails a figure left stale after a
    screening decision changed. What is drawn is the condensed flow, with each
    source's smallest reasons folded, so that is what is recomputed.
    """
    from cp_multiomics.prisma_flow import flow_table, repository_flows

    sidecar = REPO / "manuscript" / "figures" / "prisma_flow_counts.csv"
    if not sidecar.exists():
        pytest.skip("figure not drawn; run pipeline/08_figures.py --prisma-only")
    current = flow_table([flow.condensed() for flow in repository_flows(REPO)])
    drawn = pd.read_csv(sidecar, keep_default_na=False)
    # A source's pooled study units may come from results/, which a fresh clone
    # and the CI job lack; GEO's come from a tracked table and PRIDE's do not.
    # Only those rows are left out, so every other count is still checked. The
    # guard used to skip only when no source had units, which GEO's always
    # prevented: the test failed in CI while passing wherever results/ existed.
    unrecomputable = (drawn["stage"] == "units") & ~drawn["source"].isin(
        current.loc[current["stage"] == "units", "source"])
    drawn = drawn[~unrecomputable]
    pd.testing.assert_frame_equal(
        drawn.reset_index(drop=True), current.astype({"n": "int64"}).reset_index(drop=True),
        check_dtype=False, obj="prisma_flow_counts.csv (rerun 08_figures.py --prisma-only)")


def test_an_unmapped_reason_code_raises():
    """A new reason must be placed in the diagram on purpose, not lumped."""
    from cp_multiomics.prisma_flow import GEO_EXCLUDED, _group

    with pytest.raises(ValueError, match="unmapped reason"):
        _group(pd.Series(["no_pain_phenotype", "a_new_code"]), GEO_EXCLUDED, "GEO")


def test_a_flow_that_does_not_add_up_raises():
    from cp_multiomics.prisma_flow import SourceFlow

    with pytest.raises(ValueError, match="does not add up"):
        SourceFlow("X", identified=10, identified_label="studies",
                   excluded=(("reason", 3),), included=6).check()


def test_top_reasons_folds_only_when_it_saves_a_line():
    """Folding one leftover reason saves nothing and would hide its name."""
    from cp_multiomics.prisma_flow import OTHER_REASONS, SourceFlow, top_reasons

    items = (("a", 1), ("b", 9), ("c", 3), ("d", 2))
    assert top_reasons(items, 2) == (("b", 9), ("c", 3), (OTHER_REASONS, 3))
    assert top_reasons(items[:3], 2) == (("b", 9), ("c", 3), ("a", 1))
    flow = SourceFlow("X", identified=20, identified_label="studies",
                      excluded=items, included=5)
    flow.condensed(2).check()


def test_plot_repository_flow_writes_pdf_and_sidecar(tmp_path):
    from cp_multiomics.prisma_flow import OTHER_REASONS, SourceFlow
    from cp_multiomics.viz.prisma import plot_repository_flow

    flow = SourceFlow("GEO", identified=16, identified_label="series",
                      excluded=(("No contrast", 4), ("No split", 3), ("Duplicate", 2),
                                ("Unreadable", 1)),
                      awaiting=(("Pending", 1),), included=5, units=6)
    out = plot_repository_flow([flow], ["Transcriptomic: 6 units"], out_dir=tmp_path)
    assert out.exists() and out.suffix == ".pdf"
    side = pd.read_csv(tmp_path / "prisma_flow_counts.csv", keep_default_na=False)
    assert set(side["stage"]) == {"identified", "excluded", "awaiting", "included", "units"}
    # The sidecar records the fold as drawn, and the drawn counts still add up.
    excluded = side[side["stage"] == "excluded"].set_index("label")["n"]
    assert list(excluded.index) == ["No contrast", "No split", OTHER_REASONS]
    assert excluded[OTHER_REASONS] == 3


def test_a_word_wider_than_its_box_raises(tmp_path):
    """Text printing over a box edge is refused, not drawn."""
    from cp_multiomics.prisma_flow import SourceFlow
    from cp_multiomics.viz.prisma import plot_repository_flow

    flow = SourceFlow("Transcriptomicsrepositoryofrecords", identified=5,
                      identified_label="series", included=5)
    with pytest.raises(ValueError, match="wider than its box"):
        plot_repository_flow([flow], [], out_dir=tmp_path)


def test_flow_diagram_prints_at_full_width_on_one_page(tmp_path):
    """Drawn from the current records, Fig. 1 needs no scaling at print.

    It is laid out at the journal's print width with text at FONT. A screening
    change that adds a reason or wraps a label makes it taller, and past a page
    it is shrunk to fit, taking its text below FONT without any build noticing.
    """
    import re

    from cp_multiomics.prisma_flow import SUPP_DIR, repository_flows
    from cp_multiomics.viz.prisma import PAGE_HEIGHT, WIDTH, plot_repository_flow

    if not (REPO / SUPP_DIR / "Table_S2_GWAS_studies.csv").exists():
        pytest.skip("GWAS screening record absent: the code-only snapshot ships no manuscript/")
    out = plot_repository_flow(repository_flows(REPO), ["Summary line"], out_dir=tmp_path)
    box = re.search(rb"/MediaBox\s*\[\s*0\s+0\s+([\d.]+)\s+([\d.]+)\s*\]",
                    out.read_bytes())
    assert box is not None
    width, height = (float(v) / 72 for v in box.groups())
    assert width == pytest.approx(WIDTH, abs=0.01)
    assert height <= PAGE_HEIGHT, f"{height * 25.4:.0f} mm, over a page"
