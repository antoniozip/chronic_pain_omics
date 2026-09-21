#!/usr/bin/env python3
"""Emit the LaTeX body of Table 4, the cross-species sensitivity comparison.

The three arms differ only in which human studies are meta-analysed; the rodent
pools are identical across all three, which is what makes the comparison a
statement about tissue compartment rather than about pooling depth.

Generated rather than transcribed for the same reason Table 2 is: every cell is
a number that moves whenever the human corpus changes, and the arms are read
from three separate summary files that are easy to leave out of step.

Usage:
    python scripts/build_concordance_table.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]

SPECIES = [("Mus musculus", r"\emph{Mus musculus}"),
           ("Rattus norvegicus", r"\emph{Rattus norvegicus}")]


def _fmt_p(value: float) -> str:
    """Two significant figures, LaTeX scientific below 0.01."""
    if value >= 0.01:
        return f"{value:.2f}"
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}$\\times$10\\textsuperscript{{{int(exponent)}}}"


def _fmt_rho(value: float) -> str:
    text = f"{abs(value):.3f}"
    return f"$-${text}" if value < 0 else text


def arm_rows(label: str, summary: Path) -> list[str]:
    d = pd.read_csv(summary)
    out = []
    for species, printed in SPECIES:
        row = d[d["species"] == species]
        if row.empty:
            raise SystemExit(f"{summary} has no row for {species}")
        r = row.iloc[0]
        n, conc = int(r["n_ortholog_pairs"]), int(r["n_concordant"])
        # Recomputed here rather than read: the summary carries the correlation
        # p-values, not the binomial one the table reports.
        pval = stats.binomtest(conc, n, 0.5).pvalue
        out.append(f"    {label} & {printed} & {n:,} & {100 * conc / n:.1f}\\% & "
                   f"{_fmt_p(pval)} & {_fmt_rho(float(r['spearman_r']))} \\\\")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modality", default="transcriptomics")
    args = ap.parse_args()

    cross = REPO_ROOT / "results" / "cross_species" / args.modality
    sens = cross / "sensitivity"

    arms = [
        ("published_pool", sens / "published_pool" / "summary.csv"),
        ("tissue_matched", sens / "tissue_matched" / "summary.csv"),
        ("all_human", cross / "summary.csv"),
    ]
    human_n = {
        "published_pool": 5,
        "tissue_matched": 10,
        "all_human": _human_pool_size(cross),
    }
    printed = {"published_pool": "Published", "tissue_matched": "Tissue-matched",
               "all_human": "All human"}

    for key, path in arms:
        if not path.exists():
            raise SystemExit(f"missing arm summary: {path}")
        for line in arm_rows(f"{printed[key]} ({human_n[key]})", path):
            print(line)
    return 0


def _human_pool_size(cross: Path) -> int:
    """Studies in the unrestricted human pool, read from 06g's output."""
    pooled = (cross.parents[1] / "meta" / cross.name / "by_species"
              / "Homo_sapiens_pooled.csv")
    d = pd.read_csv(pooled, usecols=["study_ids"])
    return (d["study_ids"].astype(str).str.split(";").explode()
            .str.strip().nunique())


if __name__ == "__main__":
    sys.exit(main())
