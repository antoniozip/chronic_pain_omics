"""panel B must distinguish "not assayed" from "measured at zero".

`pipeline/15_gwas_tx_overlap.py` draws the pooled transcriptomic effect of the
genes genetics nominates. A gene the transcriptomic pool never measured used
to be `fillna(0.0)` and coloured as down-regulated, i.e. rendered exactly like
a genuine zero and like an effect in the opposite direction. The figure then
supported the paper's "disjoint layers" claim with an absence that looked like
a measurement.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba

from tests.conftest import load_pipeline_module


@pytest.fixture(scope="module")
def pipeline_15():
    return load_pipeline_module("15_gwas_tx_overlap.py")


def _rep(genes: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"gene": genes, "n_studies": [9] * len(genes),
                         "min_p": [1e-9] * len(genes)})


def _eff(values: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {"k": [10] * len(values), "yi_pooled": values,
         "padj_pooled": [0.01] * len(values)},
        index=pd.Index([f"G{i}" for i in range(len(values))], name="feature_id"),
    )


def test_absent_gene_is_marked_not_plotted_as_zero(pipeline_15):
    rep = _rep(["G0", "G1", "G2"])
    eff = _eff([0.8, -0.6, np.nan])

    fig, ax = plt.subplots()
    try:
        absent = pipeline_15._panel_transcriptomic(ax, rep, eff)
    finally:
        plt.close(fig)

    assert absent == 1
    assert len(ax.patches) == 3, "every gene gets a row"

    hatched = [p for p in ax.patches if p.get_hatch()]
    assert len(hatched) == 1, "the absent gene is drawn as a distinct marker"
    assert hatched[0].get_width() == 0
    assert hatched[0].get_facecolor() == to_rgba("none")

    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == ["not in pooled table (n=1)"]


def test_bar_colour_follows_the_sign_of_each_measured_gene(pipeline_15):
    rep = _rep(["G0", "G1", "G2"])
    eff = _eff([0.8, -0.6, np.nan])

    fig, ax = plt.subplots()
    try:
        pipeline_15._panel_transcriptomic(ax, rep, eff)
    finally:
        plt.close(fig)

    # Genes are drawn bottom-up: G2 (absent) at y=0, G1 (negative) at y=1,
    # G0 (positive) at y=2.
    by_row = {round(p.get_y() + p.get_height() / 2): p for p in ax.patches}
    measured = [p for p in ax.patches if not p.get_hatch()]
    assert {to_rgba(p.get_facecolor()) for p in measured} == {
        to_rgba(pipeline_15.COL_UP), to_rgba(pipeline_15.COL_DOWN)
    }
    assert to_rgba(by_row[1].get_facecolor()) == to_rgba(pipeline_15.COL_DOWN)
    assert to_rgba(by_row[2].get_facecolor()) == to_rgba(pipeline_15.COL_UP)


def test_wholly_empty_panel_fails_instead_of_drawing(pipeline_15):
    rep = _rep(["G0", "G1"])
    eff = _eff([np.nan, np.nan])

    fig, ax = plt.subplots()
    try:
        with pytest.raises(RuntimeError, match="no gene"):
            pipeline_15._panel_transcriptomic(ax, rep, eff)
    finally:
        plt.close(fig)
