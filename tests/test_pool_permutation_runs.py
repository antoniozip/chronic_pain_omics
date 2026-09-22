"""Tests for scripts/pool_permutation_runs.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pool_permutation_runs as ppr  # noqa: E402


def _run(seed: int, pcts: list[float], tested: int = 100, hits: int = 15) -> pd.DataFrame:
    """One process's table: the observed row, then its permutations."""
    rows = [{"perm": 0, "seed": None, "tested": tested, "hits": hits,
             "pct": 100 * hits / tested}]
    rows += [{"perm": i, "seed": seed * 1000 + i, "tested": 100, "hits": int(v),
              "pct": v} for i, v in enumerate(pcts, start=1)]
    return pd.DataFrame(rows)


def test_runs_pool_into_one_null_renumbered_in_seed_order():
    pooled = ppr.pool_runs({"b": _run(102, [20.0, 5.0]), "a": _run(101, [10.0, 16.0])})
    assert list(pooled["perm"]) == [0, 1, 2, 3, 4]
    assert list(pooled["seed"][1:]) == [101001, 101002, 102001, 102002]
    assert list(pooled["pct"][1:]) == [10.0, 16.0, 20.0, 5.0]


def test_p_counts_ties_as_exceedances_and_adds_one():
    pooled = ppr.pool_runs({"a": _run(101, [15.0, 20.0, 5.0, 1.0])})
    b, n, p, lo, hi = ppr.permutation_p(pooled)
    assert (b, n) == (2, 4) and p == pytest.approx(3 / 5, abs=0)
    assert 0 < lo < 0.5 < hi < 1


def test_runs_on_different_features_are_refused():
    """Their permutations are nulls for different observed statistics."""
    with pytest.raises(ValueError, match="did not test the same statistic"):
        ppr.pool_runs({"a": _run(101, [1.0]), "b": _run(102, [1.0], hits=16)})


def test_a_run_disagreeing_with_the_reference_is_refused():
    with pytest.raises(ValueError, match="did not test the same statistic"):
        ppr.pool_runs({"a": _run(101, [1.0])}, reference=(100, 16))


def test_a_repeated_seed_is_refused():
    """The same shuffle counted twice narrows the null for nothing."""
    with pytest.raises(ValueError, match="more than one permutation"):
        ppr.pool_runs({"a": _run(101, [1.0]), "b": _run(101, [2.0])})


def test_a_run_without_seeds_is_refused():
    legacy = _run(101, [1.0, 2.0]).drop(columns="seed")  # pre-seed format
    with pytest.raises(ValueError, match="without a seed"):
        ppr.pool_runs({"old": legacy, "new": _run(102, [3.0])})
