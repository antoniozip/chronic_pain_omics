"""PRISMA flow mapping tests."""

from pathlib import Path

import pandas as pd

from cp_multiomics.prisma_flow import GEO_EXCLUDED, geo_flow

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_every_excluded_candidate_reason_code_is_mapped():
    """Every reason code carried by an excluded transcriptomic candidate must
    have a GEO_EXCLUDED label, or the PRISMA flow aborts with an unmapped code.
    """
    cands = pd.read_csv(
        REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv")
    codes = set(cands.loc[cands["verdict"] == "exclude", "reason_code"].dropna())
    unmapped = sorted(codes - set(GEO_EXCLUDED))
    assert not unmapped, f"unmapped reason codes: {unmapped}"


def test_geo_flow_draws_no_zero_count_awaiting_row():
    """When the pending queue is empty, the flow must not emit awaiting rows:
    a zero in the diagram says records are still awaiting a decision."""
    flow = geo_flow(REPO_ROOT)
    assert all(n > 0 for _, n in flow.awaiting), flow.awaiting
