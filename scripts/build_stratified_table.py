#!/usr/bin/env python3
"""Emit the LaTeX body of Table 2 from the stratified meta-analysis outputs.

Table 2 showed the nine pain-model strata of the original design. The
2026-08-28 search amendment added eight more, which ``06b`` writes alongside
them, and one of those -- endometriosis (human) -- is poolable. A table that
omits it while the text says no human subtype supports meta-analytic inference
is wrong in the direction that matters, so the table now covers every stratum
``06b`` produces.

The rows are generated rather than transcribed, for the same reason Table 4 is:
Table 2's Sig column carried a stale 25 for SNI, and its Features column kept a
pre-2026-08-16 value for SNL through two green claim runs, because both were
typed by hand and read by nothing.

Ordering keeps the original nine first, then a rule, then the amendment strata.
That is the arrangement ``conf/analysis/phenotype_scope.csv`` records: the
amendment enters as its own strata and cannot move a published stratum's
numbers, and the table should show that rather than hide it.

Usage:
    python scripts/build_stratified_table.py
    python scripts/build_stratified_table.py --modality transcriptomics
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from cp_multiomics.strata import AMENDMENT, ORIGINAL, latex_label

REPO_ROOT = Path(__file__).resolve().parents[1]

# Stratum names and labels live in the package so the graphical abstract
# prints the same words this table does; see cp_multiomics/strata.py.

logger = logging.getLogger(__name__)


def _fmt_padj(value: float) -> str:
    """Format a padj the way the table's existing cells are written."""
    if pd.isna(value):
        return "---"
    if value == 0:
        # A bound cannot be mutation-tested: anything satisfies "< x". Say so
        # loudly rather than printing "$<$1e-300" as the old CIPN cell did.
        return r"underflow"
    # Python's default two-digit exponent is the convention already in the
    # table (1.1e-20, 1.2e-147), and the claim patterns match on it.
    return f"{value:.1e}"


def stratum_row(path: Path, label: str) -> str:
    """One LaTeX row: Model, Studies, Features, Pooled, Sig, Top gene (padj)."""
    d = pd.read_csv(path, usecols=["feature_id", "k", "study_ids", "pval",
                                   "padj"])
    studies = (d["study_ids"].astype(str).str.split(";").explode()
               .str.strip().nunique())
    family = d[d["padj"].notna()]
    n_sig = int((family["padj"] < 0.05).sum())

    if len(family):
        # Rank on padj, then raw p: BH assigns tied adjusted values to adjacent
        # ranks, so padj alone leaves the top cell arbitrary. SNI's top two
        # both sit at 5.0e-05 and are separated only here.
        top = family.sort_values(["padj", "pval"]).iloc[0]
        gene = f"{top['feature_id']} ({_fmt_padj(float(top['padj']))})"
    else:
        gene = "---"

    return (f"    {label} & {studies} & {len(d):,} & {len(family):,} & "
            f"{n_sig:,} & {gene} \\\\")


def q_row(path: Path) -> str:
    """The Between-model Q row: features tested, and how many differ by model.

    Hand-written until 2026-08-31, and stale by two regenerations when it was
    caught -- it read 33,853 tested and 3,685 significant while the sentence
    quoting the same test said 62,792 and 13,192. Neither number was claimed,
    the Discussion restated the 3,685, and nothing compared the two.
    """
    d = pd.read_csv(path, usecols=["feature_id", "QM_pval", "padj_QM"])
    n_sig = int((d["padj_QM"] < 0.05).sum())
    top = d.sort_values(["padj_QM", "QM_pval"]).iloc[0]["feature_id"]
    return (f"    Between-model Q & --- & {len(d):,} & --- & "
            f"\\textbf{{{n_sig:,}}} & {top} \\\\")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modality", default="transcriptomics")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    strat = (REPO_ROOT / "results" / "meta" / args.modality / "stratified")
    if not strat.exists():
        raise SystemExit(f"no stratified results at {strat}")

    known = {stem for stem, _ in ORIGINAL + AMENDMENT}
    on_disk = {p.name.replace("_pooled.csv", "")
               for p in strat.glob("*_pooled.csv")}
    missing = sorted(on_disk - known)
    if missing:
        # A stratum 06b writes but this script does not know about would be
        # silently absent from the paper -- the failure being repaired here.
        raise SystemExit(
            f"06b produced strata this table does not list: {missing}. "
            "Add them to ORIGINAL or AMENDMENT before regenerating."
        )

    for group in (ORIGINAL, AMENDMENT):
        for stem, _ in group:
            path = strat / f"{stem}_pooled.csv"
            if not path.exists():
                logger.warning("no table for %s; skipping", stem)
                continue
            print(stratum_row(path, latex_label(stem)))
        if group is ORIGINAL:
            print(r"    \midrule")

    q_path = strat / "between_model_Q.csv"
    if q_path.exists():
        print(r"    \midrule")
        print(q_row(q_path))
    else:
        logger.warning("no between_model_Q.csv; omitting the Q row")
    return 0


if __name__ == "__main__":
    sys.exit(main())
