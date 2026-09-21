#!/usr/bin/env python
"""Re-key PXD013362 from peptide sequences to gene symbols.

PXD013362 is a neuropeptidomics study: its features are peptide sequences,
while PXD054342 and PXD055816 are keyed on gene symbols. Two feature spaces
means the three studies never meet, so `meta.min_studies = 3` yields an empty
correction family however good each study is. This maps one space onto the
other.

The mapping itself lives in `conf/analysis/peptide_gene_map.csv` rather than in
code, because assigning a peptide to its precursor is a biological judgement
that should be reviewable and citable. Rows with an empty `peptide` are
documented targets, not assignments, and are ignored here — they record which
precursors the study reports (CALCA, ADCYAP1, VIP, SCG2, SCG3, PENK, TRH,
PCSK1N) so the list is in the repository while the per-peptide assignments are
filled in from the submitter's annotation or a sequence search.

**One peptide represents a gene.** Several peptides from one precursor,
measured in the same animals, are correlated; pooling them as independent
observations would reintroduce exactly the pseudo-replication that collapsing
the brain regions removed. The peptide with the smallest p-value is kept, which
is the rule `05b`'s `aggregate_to_gene` already applies for transcriptomics.

Input:  results/per_study/proteomics/PXD013362/PXD013362_effects_collapsed.csv
Output: results/per_study/proteomics/PXD013362/PXD013362_effects_gene.csv

Usage:
    python pipeline/05i_peptide_to_gene.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
MAP = REPO_ROOT / "conf" / "analysis" / "peptide_gene_map.csv"
STUDY_DIR = REPO_ROOT / "results" / "per_study" / "proteomics" / "PXD013362"
IN = STUDY_DIR / "PXD013362_effects_collapsed.csv"
OUT = STUDY_DIR / "PXD013362_effects_gene.csv"


def load_map() -> tuple[dict[str, str], list[str]]:
    """Return peptide -> gene assignments, and the documented target genes."""
    if not MAP.exists():
        raise SystemExit(f"no mapping table at {MAP}")
    m = pd.read_csv(MAP, dtype=str).fillna("")
    m["peptide"] = m["peptide"].str.strip()
    m["gene_symbol"] = m["gene_symbol"].str.strip()

    assigned = m[m["peptide"] != ""]
    targets = sorted(set(m.loc[m["gene_symbol"] != "", "gene_symbol"]))
    dupes = assigned["peptide"][assigned["peptide"].duplicated()].tolist()
    if dupes:
        raise SystemExit(f"peptide assigned to more than one gene: {dupes[:5]}")
    return dict(zip(assigned["peptide"], assigned["gene_symbol"])), targets


def main() -> int:
    mapping, targets = load_map()
    print(f"mapping table: {len(mapping)} peptide assignments, "
          f"{len(targets)} target genes ({', '.join(targets)})")

    if not IN.exists():
        raise SystemExit(
            f"no collapsed table at {IN}\n"
            "Run pipeline/05h_collapse_pxd013362.R first: peptides must be "
            "collapsed across brain regions before they are re-keyed, or the "
            "region replicates would be carried into the gene-level table."
        )

    df = pd.read_csv(IN)
    df["feature_id"] = df["feature_id"].astype(str)

    if not mapping:
        print(
            "\nNo peptide assignments yet, so nothing was written.\n"
            "Fill the `peptide` column of conf/analysis/peptide_gene_map.csv,\n"
            "one row per peptide sequence, taking the assignments from the\n"
            "submitter's annotation in peptide_quantitation.xlsx or from a\n"
            "search of each sequence against the mouse precursors listed there.\n"
            f"The table to map against has {df['feature_id'].nunique()} distinct "
            "peptides."
        )
        return 0

    df["gene_symbol"] = df["feature_id"].map(mapping)
    mapped = df[df["gene_symbol"].notna()].copy()
    print(f"\n{mapped['feature_id'].nunique()} of {df['feature_id'].nunique()} "
          f"peptides mapped, covering {mapped['gene_symbol'].nunique()} genes")
    if mapped.empty:
        print("no peptide in the table matched an assignment; nothing written")
        return 0

    # One peptide per gene per study unit: the most significant. See module
    # docstring for why pooling peptides of a precursor is not an option.
    mapped = mapped.sort_values("pval")
    best = mapped.drop_duplicates(subset=["gene_symbol", "study_id"], keep="first")
    best = best.assign(
        original_feature_id=best["feature_id"], feature_id=best["gene_symbol"]
    ).drop(columns=["gene_symbol"])

    best.to_csv(OUT, index=False)
    per_unit = best.groupby("study_id")["feature_id"].nunique().to_dict()
    print(f"-> {OUT.relative_to(REPO_ROOT)}  {len(best)} rows")
    print(f"   genes per study unit: {per_unit}")
    print("\nThese genes reach k=3 only if PXD054342 and PXD055816 also measured "
          "them; CALCA is the one worth checking first.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
