#!/usr/bin/env python3
"""Emit the LaTeX body of Table 1, the leading well-replicated pooled genes.

Ranking the 773 significant features by adjusted p-value alone puts the
shallowest ones first: the ten that led the former forest-plot figure were
pooled from a median of five study units, most of them non-coding or predicted
loci. Table 1 instead ranks within the features pooled from at least as many
study units as the significant set's median, which is computed here, not
fixed, so the rule follows the analysis when it moves.

Every cell except the role is read from ``pooled_effects.csv``.
`tests/test_manuscript_claims.py::test_top_genes_table_matches_the_pool`
regenerates the rows and compares them with the manuscript, so a
regeneration of the pool fails the suite until this script is re-run and
its output pasted in.

The roles are short functional annotations (NCBI Gene / UniProt summaries),
not results of this analysis. A gene that enters the table without one is an
error, not a blank cell.

Usage:
    python scripts/build_top_genes_table.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
POOLED = REPO_ROOT / "results" / "meta" / "transcriptomics" / "pooled_effects.csv"
N_ROWS = 10

ROLES: dict[str, str] = {
    "Atf3": "Injury-induced transcription factor; marker of axotomised "
            "sensory neurons",
    "Tspo": "Mitochondrial translocator protein; marker of glial activation",
    "Clic1": "Intracellular chloride channel; microglial activation",
    "P2ry6": "Microglial UDP receptor; drives phagocytosis",
    "Timp1": "Inhibitor of matrix metalloproteinases; tissue remodelling "
             "after injury",
    "Icam1": "Endothelial adhesion molecule; leukocyte recruitment",
    "Fcgr2b": "Inhibitory Fc$\\gamma$ receptor; immune-complex signalling",
    "C3": "Central complement component",
    "Cd44": "Hyaluronan receptor; immune-cell adhesion and migration",
    "Sbno2": "IL-10/STAT3-induced regulator; restrains macrophage "
             "inflammation",
}


def depth_threshold(pooled: pd.DataFrame) -> int:
    """Median study-unit count of the significant set (padj < 0.05)."""
    sig = pooled[pooled["padj_pooled"] < 0.05]
    return int(sig["k"].median())


def top_genes(pooled: pd.DataFrame, n: int = N_ROWS) -> pd.DataFrame:
    """The ``n`` most significant features pooled from >= the median depth.

    Ranked on padj and then raw p, as the other generated tables are, so
    tied adjusted values keep a defined order.
    """
    k_min = depth_threshold(pooled)
    sig = pooled[(pooled["padj_pooled"] < 0.05) & (pooled["k"] >= k_min)]
    return sig.sort_values(["padj_pooled", "pval_pooled"]).head(n)


def n_series(study_ids: str) -> int:
    """Distinct GEO series behind the study units (a split unit is one series)."""
    return len({re.match(r"GSE\d+", s.strip()).group(0)  # type: ignore[union-attr]
                for s in study_ids.split(";")})


def _fmt_signed(value: float) -> str:
    text = f"{abs(value):.2f}"
    return f"$-${text}" if value < 0 else text


def _fmt_p(value: float) -> str:
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}$\\times$10\\textsuperscript{{{int(exponent)}}}"


def _fmt_i2(value: float) -> str:
    return "$>$99.9" if value >= 99.95 else f"{value:.1f}"


def table_rows(pooled: pd.DataFrame) -> list[str]:
    rows = []
    for _, r in top_genes(pooled).iterrows():
        gene = str(r["feature_id"])
        if gene not in ROLES:
            raise SystemExit(f"no role written for {gene}; add it to ROLES")
        rows.append(
            f"    \\gene{{{gene}}} & {n_series(str(r['study_ids']))} "
            f"({int(r['k'])}) & {_fmt_signed(r['yi_pooled'])} & "
            f"{_fmt_signed(r['ci_lb'])} to {_fmt_signed(r['ci_ub'])} & "
            f"{_fmt_i2(r['I2'])} & {_fmt_p(r['padj_pooled'])} & "
            f"{ROLES[gene]} \\\\")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pooled", type=Path, default=POOLED)
    args = ap.parse_args()
    pooled = pd.read_csv(args.pooled)
    for line in table_rows(pooled):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
