"""Tests for viz package and pipeline/08_figures.py."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")  # non-interactive backend for tests

from cp_multiomics.viz.concordance_plot import plot_concordance_scatter
from cp_multiomics.viz.heatmap import _cluster_order, plot_cross_modal_heatmap
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
# Cross-modal heatmap
# ---------------------------------------------------------------------------

def test_plot_cross_modal_heatmap_creates_pdf(tmp_path):
    data = {
        "transcriptomics": _pooled_df(seed=1),
        "proteomics":      _pooled_df(seed=2),
    }
    out = plot_cross_modal_heatmap(data, out_dir=tmp_path, top_n=10, padj_threshold=0.05)
    assert out.exists()


def test_heatmap_no_significant_features_returns_path(tmp_path):
    data = {
        "transcriptomics": _pooled_df().assign(padj_pooled=1.0),
        "proteomics":      _pooled_df().assign(padj_pooled=1.0),
    }
    out = plot_cross_modal_heatmap(data, out_dir=tmp_path)
    # Returns path even when empty (may be a .txt sentinel)
    assert out is not None


def test_cluster_order_returns_correct_length():
    data = np.random.rand(5, 3)
    order = _cluster_order(data)
    assert len(order) == 5
    assert sorted(order) == list(range(5))


def test_cluster_order_single_row():
    data = np.array([[1.0, 2.0, 3.0]])
    order = _cluster_order(data)
    assert order == [0]


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
    screening decision changed.
    """
    from cp_multiomics.prisma_flow import flow_table, repository_flows

    sidecar = REPO / "manuscript" / "figures" / "prisma_flow_counts.csv"
    if not sidecar.exists():
        pytest.skip("figure not drawn; run pipeline/08_figures.py --prisma-only")
    current = flow_table(repository_flows(REPO))
    drawn = pd.read_csv(sidecar, keep_default_na=False)
    if "units" in set(drawn["stage"]) - set(current["stage"]):
        pytest.skip("results/ absent, so the pooled study units cannot be recomputed")
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


def test_plot_repository_flow_writes_pdf_and_sidecar(tmp_path):
    from cp_multiomics.prisma_flow import SourceFlow
    from cp_multiomics.viz.prisma import plot_repository_flow

    flow = SourceFlow("GEO", identified=10, identified_label="series",
                      excluded=(("No contrast", 4),), awaiting=(("Pending", 1),),
                      included=5, units=6)
    out = plot_repository_flow([flow], ["Transcriptomic: 6 units"], out_dir=tmp_path)
    assert out.exists() and out.suffix == ".pdf"
    side = pd.read_csv(tmp_path / "prisma_flow_counts.csv", keep_default_na=False)
    assert set(side["stage"]) == {"identified", "excluded", "awaiting", "included", "units"}
