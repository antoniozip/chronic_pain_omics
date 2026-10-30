"""Per-study BH adjustment must survive a non-finite p-value.

`multipletests(..., method="fdr_bh")` sorts NaN last and then
`np.minimum.accumulate` propagates it across the whole array, so ONE NaN
p-value makes every `padj` in that study NaN while `reject` still reports the
finite ones significant. Downstream `(padj < 0.05).sum()` is then 0: the study
silently contributes zero significant features.

Both per-study DA stages now adjust over the non-NA p-values and leave the
NA-p rows NA.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

PIPELINE = Path(__file__).parent.parent / "pipeline"


def _load(filename: str):
    name = "_test_" + filename.replace(".", "_").replace("-", "_")
    spec = importlib.util.spec_from_file_location(name, PIPELINE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def proteomics():
    return _load("05g_proteomics_pride_da.py")


@pytest.fixture(scope="module")
def metabolomics():
    return _load("05j_metabolomics_da.py")


P = np.array([0.01, 0.02, 0.03, np.nan, 0.5])
FINITE = [0, 1, 2, 4]


def test_05g_bh_adjust_keeps_the_family_when_one_pval_is_nan(proteomics):
    got = proteomics.bh_adjust(P)
    expected = multipletests(P[FINITE], method="fdr_bh")[1]
    np.testing.assert_allclose(got[FINITE], expected)
    assert np.isnan(got[3])
    assert np.isfinite(got[FINITE]).all()


def test_05g_bh_adjust_matches_the_clean_family_exactly(proteomics):
    got = proteomics.bh_adjust(P)
    clean = proteomics.bh_adjust(P[FINITE])
    np.testing.assert_allclose(got[FINITE], clean)


def test_05j_bh_adjust_keeps_the_family_when_one_pval_is_nan(metabolomics):
    got = metabolomics.bh_adjust(P)
    expected = multipletests(P[FINITE], method="fdr_bh")[1]
    np.testing.assert_allclose(got[FINITE], expected)
    assert np.isnan(got[3])
    assert np.isfinite(got[FINITE]).all()


def test_05j_bh_adjust_matches_the_clean_family_exactly(metabolomics):
    got = metabolomics.bh_adjust(P)
    clean = metabolomics.bh_adjust(P[FINITE])
    np.testing.assert_allclose(got[FINITE], clean)


def test_05g_differential_adjusts_through_the_helper(proteomics):
    import inspect

    assert "bh_adjust(" in inspect.getsource(proteomics.differential)


def test_05j_finalise_adjusts_through_the_helper(metabolomics):
    import inspect

    assert "bh_adjust(" in inspect.getsource(metabolomics.finalise)
