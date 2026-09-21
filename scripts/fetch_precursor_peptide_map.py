#!/usr/bin/env python
"""Assign PXD013362 peptides to precursor genes by substring match.

PXD013362 is keyed on peptide sequences; the other two proteomics studies are
keyed on gene symbols. `pipeline/05i_peptide_to_gene.py` bridges the two using
`conf/analysis/peptide_gene_map.csv`, which this script fills.

The assignment rule is deliberately narrow: a peptide belongs to a precursor if
its bare amino-acid sequence is a **substring** of that precursor's mouse
UniProt sequence. Peptidomics reports endogenous cleavage products, so a peptide
of a neuropeptide precursor appears verbatim in it — no enzymatic-digest
reasoning is needed. Anything matching two of the eight precursors is left
unassigned rather than guessed at.

Two things the peptide strings need before matching:

- **Modification tokens are stripped.** The study writes them inline in two
  notations, `A(+42.01)ADISQW` and `A[+42]ADISQW`; both are removed.
- **Short peptides are rejected unless UniProt names them.** A 6-mer occurs in a
  617-residue precursor by chance often enough to matter (~0.3% per precursor at
  uniform composition, and far more for the low-complexity stretches these
  precursors are full of), so `MIN_LEN` is the length below which a substring
  match is not evidence on its own. A short sequence that exactly equals an
  annotated `PEPTIDE` feature of the precursor is kept anyway: that is
  independent curated evidence, not a coincidence of length. This is what
  rescues Met-enkephalin (`YGGFM`) and Leu-enkephalin (`YGGFL`) — the two
  best-characterised products in the whole set, and both 5-mers. It does not
  rescue `GGFMRF`, which is a truncation of Met-enkephalin-Arg-Phe rather than
  an annotated product itself.

Matching is restricted to eight precursors, so a match says "this sequence is
present in gene X", not "this peptide is unique to gene X in the mouse
proteome". That is acceptable here because the eight are the precursors the
study itself reports, and because the alternative — an unmapped peptide — costs
the feature entirely.

Usage:
    python scripts/fetch_precursor_peptide_map.py            # dry run
    python scripts/fetch_precursor_peptide_map.py --apply    # write the map
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
MAP = REPO_ROOT / "conf" / "analysis" / "peptide_gene_map.csv"
COLLAPSED = (
    REPO_ROOT
    / "results"
    / "per_study"
    / "proteomics"
    / "PXD013362"
    / "PXD013362_effects_collapsed.csv"
)

UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
UNIPROT_FEATURES = "https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ft_peptide"

# One reviewed mouse (organism 10090) entry per target gene. Resolved by
# `gene:<symbol> AND organism_id:10090 AND reviewed:true` on 2026-08-15.
#
# Keyed on the **MGI mouse symbol**, not the uppercase human one. PXD013362 is a
# mouse study, and the two studies it has to pool with (PXD054342, PXD055816)
# write `Calca`/`Penk`/`Scg2`. Every join downstream is a case-sensitive string
# match, so emitting `CALCA` here would leave all eight genes at k=1 without
# erroring anywhere — the silent-fragmentation failure `R/symbol_utils.R` exists
# to prevent. Human and rodent symbols are separate features throughout this
# pipeline, joined only through the ortholog map.
PRECURSORS: dict[str, tuple[str, str]] = {
    "Calca": ("Q99JA0", "calcitonin gene-related peptide (CGRP)"),
    "Adcyap1": ("O70176", "pituitary adenylate cyclase-activating polypeptide (PACAP)"),
    "Vip": ("P32648", "vasoactive intestinal peptide"),
    "Scg2": ("Q03517", "secretogranin II"),
    "Scg3": ("P47867", "secretogranin III"),
    "Penk": ("P22005", "proenkephalin"),
    "Trh": ("Q62361", "thyrotropin-releasing hormone"),
    "Pcsk1n": ("Q9QXV0", "proSAAS / proprotein convertase subtilisin-kexin type 1 inhibitor"),
}

# The gene-keyed studies PXD013362 must pool with. Checked, not assumed: a
# casing or nomenclature drift between the map and these tables produces k=1
# everywhere rather than an error.
PARTNER_STUDIES = ("PXD054342", "PXD055816")

# Calcitonin is the other splice product of Calca. It shares its N-terminus with
# CGRP1 (Q99JA0) and diverges after it, so it can only add peptides, never move
# one to a different gene. Checked as a diagnostic; see `--include-calcitonin`.
CALCA_CALCITONIN = "P70160"

MIN_LEN = 7

MOD_TOKEN = re.compile(r"\([^)]*\)|\[[^]]*\]")


def bare_sequence(peptide: str) -> str:
    """Strip inline modification tokens, leaving amino acids only."""
    return MOD_TOKEN.sub("", peptide).strip().upper()


def fetch_sequence(accession: str) -> str:
    """Return the amino-acid sequence of a UniProt entry."""
    resp = requests.get(UNIPROT.format(acc=accession), timeout=60)
    resp.raise_for_status()
    lines = resp.text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SystemExit(f"{accession}: not a FASTA response")
    return "".join(line.strip() for line in lines[1:])


def fetch_named_peptides(accession: str, sequence: str) -> dict[str, str]:
    """Return {subsequence: description} for each annotated PEPTIDE feature.

    UniProt gives cleavage products as 1-based inclusive coordinates into the
    precursor, so the sequence has to be sliced out rather than read directly.
    """
    resp = requests.get(UNIPROT_FEATURES.format(acc=accession), timeout=60)
    resp.raise_for_status()
    named: dict[str, str] = {}
    for feat in resp.json().get("features", []):
        if feat.get("type") != "Peptide":
            continue
        loc = feat["location"]
        start, end = loc["start"].get("value"), loc["end"].get("value")
        if start is None or end is None:
            continue
        named[sequence[start - 1:end]] = feat.get("description") or "unnamed peptide"
    return named


def check_partner_studies(genes: list[str]) -> None:
    """Report how many partner studies carry each assigned symbol, verbatim.

    A gene reaches the `meta.min_studies = 3` correction family only if the
    other two studies wrote it the same way, so this is the number that decides
    whether the mapping bought anything.
    """
    per_study = REPO_ROOT / "results" / "per_study" / "proteomics"
    present: dict[str, list[str]] = {}
    for acc in PARTNER_STUDIES:
        path = per_study / acc / f"{acc}_effects.csv"
        if not path.exists():
            print(f"\n[warn] {path.relative_to(REPO_ROOT)} absent; cannot verify overlap")
            continue
        features = set(pd.read_csv(path, usecols=["feature_id"])["feature_id"].astype(str))
        present[acc] = [g for g in genes if g in features]

    if not present:
        return
    print("\nOverlap with the gene-keyed studies (exact string match):")
    for gene in genes:
        carriers = [acc for acc, hits in present.items() if gene in hits]
        # PXD013362 contributes two independent units (Migraine, OIH).
        k = len(carriers) + 2
        print(f"  {gene:8s} k={k}  "
              f"{'+'.join(carriers) if carriers else 'neither partner study'}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="write conf/analysis/peptide_gene_map.csv")
    ap.add_argument("--include-calcitonin", action="store_true",
                    help="also match against the calcitonin splice product of Calca")
    args = ap.parse_args()

    if not COLLAPSED.exists():
        raise SystemExit(
            f"no collapsed table at {COLLAPSED}\n"
            "Run pipeline/05h_collapse_pxd013362.R first."
        )

    targets = {gene: (acc, prec) for gene, (acc, prec) in PRECURSORS.items()}
    sequences: dict[str, list[tuple[str, str]]] = {}
    named: dict[str, dict[str, str]] = {}
    for gene, (acc, _prec) in targets.items():
        seq = fetch_sequence(acc)
        sequences[gene] = [(acc, seq)]
        named[acc] = fetch_named_peptides(acc, seq)
        print(f"{gene:8s} {acc}  {len(seq)} aa, "
              f"{len(named[acc])} annotated cleavage products")
    if args.include_calcitonin:
        seq = fetch_sequence(CALCA_CALCITONIN)
        sequences["CALCA"].append((CALCA_CALCITONIN, seq))
        named[CALCA_CALCITONIN] = fetch_named_peptides(CALCA_CALCITONIN, seq)
        print(f"{'CALCA':8s} {CALCA_CALCITONIN}  {len(seq)} aa  (calcitonin isoform)")

    df = pd.read_csv(COLLAPSED)
    peptides = sorted(df["feature_id"].astype(str).unique())
    print(f"\n{len(peptides)} distinct peptides in {COLLAPSED.name}")

    rows: list[dict[str, str]] = []
    ambiguous: list[tuple[str, list[str]]] = []
    too_short: list[str] = []
    rescued: list[str] = []
    for peptide in peptides:
        bare = bare_sequence(peptide)
        hits = [
            (gene, acc)
            for gene, entries in sequences.items()
            for acc, seq in entries
            if bare in seq
        ]
        genes = sorted({gene for gene, _ in hits})
        if not genes:
            continue
        if len(genes) > 1:
            ambiguous.append((peptide, genes))
            continue
        gene = genes[0]
        acc = next(a for g, a in hits if g == gene)

        # A short match is a coincidence unless UniProt names it. See module
        # docstring: this is what keeps Met- and Leu-enkephalin.
        annotation = named.get(acc, {}).get(bare)
        if len(bare) < MIN_LEN and annotation is None:
            too_short.append(f"{peptide} ({len(bare)} aa -> {gene})")
            continue
        evidence = (f"named_cleavage_product_of_uniprot:{acc}:{annotation}"
                    if annotation is not None
                    else f"substring_of_uniprot:{acc}")
        rows.append({
            "peptide": peptide,
            "gene_symbol": gene,
            "precursor": targets[gene][1],
            "evidence": evidence,
        })
        if annotation is not None:
            rescued.append(f"{peptide} ({len(bare)} aa) = {gene} {annotation}")

    assigned = pd.DataFrame(rows)
    print(f"\n{len(assigned)} peptides assigned to "
          f"{assigned['gene_symbol'].nunique() if not assigned.empty else 0} genes")
    if not assigned.empty:
        for gene, n in assigned["gene_symbol"].value_counts().sort_index().items():
            print(f"  {gene:8s} {n}")
    if rescued:
        print(f"\n{len(rescued)} assignments carry a UniProt cleavage annotation:")
        for line in rescued:
            print(f"  {line}")
    if too_short:
        print(f"\n{len(too_short)} matches rejected as shorter than {MIN_LEN} aa "
              "with no UniProt cleavage annotation:")
        for line in too_short[:10]:
            print(f"  {line}")
    if ambiguous:
        print(f"\n{len(ambiguous)} peptides matched more than one precursor, "
              "left unassigned:")
        for peptide, genes in ambiguous[:10]:
            print(f"  {peptide} -> {','.join(genes)}")

    # Documented targets keep their rows with an empty `peptide`, so the eight
    # precursors the study reports stay on the record whether or not any of
    # their peptides survived matching. 05i ignores those rows.
    documented = pd.DataFrame([
        {"peptide": "", "gene_symbol": gene, "precursor": prec,
         "evidence": "target_gene_awaiting_peptide_assignment"}
        for gene, (_acc, prec) in targets.items()
        if assigned.empty or gene not in set(assigned["gene_symbol"])
    ])
    out = pd.concat([documented, assigned], ignore_index=True)

    check_partner_studies(sorted(set(assigned["gene_symbol"])) if not assigned.empty else [])

    if not args.apply:
        print(f"\nDry run. {len(out)} rows would be written to "
              f"{MAP.relative_to(REPO_ROOT)}; pass --apply to write.")
        return 0

    out.to_csv(MAP, index=False)
    print(f"\n-> {MAP.relative_to(REPO_ROOT)}  {len(out)} rows "
          f"({len(assigned)} assignments, {len(documented)} unassigned targets)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
