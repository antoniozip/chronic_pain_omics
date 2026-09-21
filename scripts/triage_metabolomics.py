#!/usr/bin/env python
"""Open every included MetaboLights study and record what it actually ships.

Screening (`literature/prisma/metabolomics_manual_review.csv`) asked whether a
study *sounds* like a pain-versus-control design. This asks the only question
that decides whether it can enter a meta-analysis: does it deposit a
quantification table with per-sample values, identified metabolites, and group
labels that encode the contrast?

MetaboLights studies are ISA-Tab. Three files matter:

- ``m_*_maf.tsv``  the metabolite assignment file: one row per metabolite, the
  fixed MAF columns, then one abundance column per sample. This is the
  quantification table; without it there is nothing to pool.
- ``s_*.txt``      the sample table, whose ``Factor Value[...]`` columns carry
  the experimental groups. A MAF with no recoverable grouping is the
  "usable with work" case from the proteomic triage.
- ``a_*.txt``      the assay table, linking samples to raw files. Counted only
  to report how many assays a study spreads itself over.

**Identified metabolites are the pooling currency**, and identification comes
at two grades that must not be conflated. ``database_identifier`` holds a ChEBI
accession and is the join key the pipeline wants. ``metabolite_identification``
holds a name, which is weaker — names need resolving to ChEBI before two
studies can be joined — but it is not nothing.

Judging on ``database_identifier`` alone is wrong, and wrong in the expensive
direction: MTBLS14059 deposits an empty ChEBI column beside names like
"L-Tyrosine" and "Tryptophan". Discarding it as unidentified would have thrown
away a real quantification table over a missing accession that a name lookup
recovers. So a row counts as identified if it carries either, and the two are
reported separately.

The opposite error is just as easy. A name can be a placeholder rather than an
identification — ``Cluster_070881`` in MTBLS12704, ``metabolite1`` in
MTBLS14319 — which is an unannotated feature wearing a name column.
``UNNAMED`` matches those, and a study whose names are *all* placeholders is
genuinely unusable. This is the metabolomic form of the peptide-versus-gene
identifier problem that blocked the proteomic arm.

Usage:
    python scripts/triage_metabolomics.py
    python scripts/triage_metabolomics.py --accessions MTBLS13869,MTBLS9662
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import re
import sys
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
REVIEW = REPO_ROOT / "literature" / "prisma" / "metabolomics_manual_review.csv"
OUT = REPO_ROOT / "literature" / "prisma" / "metabolomics_file_triage.csv"

FILES_API = "https://www.ebi.ac.uk/metabolights/ws/studies/{acc}/files"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{acc}/{name}"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# Fixed MAF columns across MS and NMR variants; anything else in the header is
# taken to be a per-sample abundance column.
MAF_FIXED = {
    "database_identifier", "chemical_formula", "smiles", "inchi",
    "metabolite_identification", "mass_to_charge", "fragmentation",
    "modifications", "charge", "retention_time", "taxid", "species",
    "database", "database_version", "reliability", "uri", "search_engine",
    "search_engine_score", "smallmolecule_abundance_sub",
    "smallmolecule_abundance_stdev_sub", "smallmolecule_abundance_std_error_sub",
    "chemical_shift", "multiplicity",
}

# Values that occupy `metabolite_identification` without identifying anything:
# feature cluster ids, bare indices, and placeholder names.
UNNAMED = re.compile(
    r"^(cluster[_\- ]?\d+|unknown.*|unidentified.*|peak[_\- ]?\d+|m\d+t\d+"
    r"|feature[_\- ]?\d+|metabolite[_\- ]?\d+|compound[_\- ]?\d+|\d+(\.\d+)?)$",
    re.I,
)

# A factor level naming the comparator arm. Two or more levels is not by itself
# a case-versus-control contrast: MTBLS9662 has Mild/Moderate/Severe, 30 each,
# which is a severity gradient among cases with no unaffected group. Pooling
# that as if it were case-versus-control would compare severe disease against
# mild disease and call the result a disease effect.
CONTROL_LEVEL = re.compile(
    r"\b(control|normal|healthy|sham|naive|na[iy]ve|vehicle|baseline|untreated"
    r"|wild[_\- ]?type|wt|non[_\- ]?(pain|disease|oa)|placebo)\b",
    re.I,
)

# Levels that are apparatus, not biology, and must not be read as an arm.
NON_BIOLOGICAL_LEVEL = re.compile(r"^(qc|quality control|blank|pool(ed)?.*)$", re.I)


def get(url: str) -> requests.Response:
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    return r


def list_files(acc: str) -> list[dict]:
    return get(FILES_API.format(acc=acc)).json().get("study", []) or []


def read_tsv(acc: str, name: str) -> tuple[list[str], list[list[str]]]:
    text = get(FTP.format(acc=acc, name=name)).content.decode("utf-8", "replace")
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if not rows:
        return [], []
    return [c.strip() for c in rows[0]], rows[1:]


def inspect_maf(acc: str, name: str) -> dict:
    header, rows = read_tsv(acc, name)
    sample_cols = [c for c in header if c and c not in MAF_FIXED]

    # A declared sample column is not a measured one. MTBLS13513 declares 91
    # columns beside 1,225 named metabolites and fills none of them; so do
    # MTBLS14059 and MTBLS14295. Counting columns rather than values called all
    # three usable, and the DA then produced nothing at all from them. What
    # makes a study poolable is populated columns, so that is what is counted.
    idx = {c: i for i, c in enumerate(header)}
    populated = []
    for c in sample_cols:
        i = idx[c]
        for r in rows:
            if len(r) > i and r[i].strip():
                try:
                    float(r[i])
                except ValueError:
                    continue
                populated.append(c)
                break

    def col(field: str) -> int:
        return header.index(field) if field in header else -1

    di, mi = col("database_identifier"), col("metabolite_identification")
    with_db = with_name = 0
    for r in rows:
        if di >= 0 and len(r) > di and r[di].strip():
            with_db += 1
        if mi >= 0 and len(r) > mi:
            v = r[mi].strip()
            if v and not UNNAMED.match(v):
                with_name += 1
    return {
        "maf": name,
        "metabolites": len(rows),
        # A ChEBI accession is the join key; a real name is recoverable to one.
        # Either counts as identified, and both are reported.
        "identified": max(with_db, with_name),
        "with_db_id": with_db,
        "with_name": with_name,
        "sample_columns": len(populated),
        "declared_columns": len(sample_cols),
    }


def inspect_samples(acc: str, name: str) -> dict:
    header, rows = read_tsv(acc, name)
    factors: dict[str, dict[str, int]] = {}
    for i, col in enumerate(header):
        if not col.startswith("Factor Value["):
            continue
        label = col[len("Factor Value["):].rstrip("]")
        levels: dict[str, int] = {}
        for r in rows:
            if len(r) > i and r[i].strip():
                levels[r[i].strip()] = levels.get(r[i].strip(), 0) + 1
        if levels:
            factors[label] = levels
    return {"n_samples": len(rows), "factors": factors}


def verdict(maf: dict | None, samples: dict) -> tuple[str, str]:
    """Classify a study, returning (verdict, reason_code)."""
    if maf is None:
        return "unusable", "no_quantification_table"
    if maf["sample_columns"] == 0:
        if maf.get("declared_columns"):
            # Columns are declared and every one is empty: identifications
            # deposited without the abundances they identify.
            return "unusable", "maf_columns_declared_but_empty"
        return "unusable", "maf_without_per_sample_values"
    if maf["identified"] == 0:
        return "unusable", "unannotated_features_only"
    if not samples["factors"]:
        return "usable_with_work", "no_factor_value_groups"
    # A contrast needs a factor with two or more *biological* levels, one of
    # which is a comparator arm. QC and blanks are apparatus, and a binary
    # 0/1 factor encodes its own control.
    has_contrast = False
    for levels in samples["factors"].values():
        bio = [lv for lv in levels if not NON_BIOLOGICAL_LEVEL.match(lv)]
        if len(bio) < 2:
            continue
        if set(bio) == {"0", "1"} or any(CONTROL_LEVEL.search(lv) for lv in bio):
            has_contrast = True
            break
    if not has_contrast:
        return "usable_with_work", "no_control_arm_in_factors"
    if maf["with_db_id"] == 0:
        # Names without ChEBI: joinable across studies only after a name lookup,
        # so it is real work, not a formality — synonyms and stereochemistry
        # make name matching lossy in exactly the way accessions avoid.
        return "usable_with_work", "named_metabolites_without_database_identifier"
    return "usable", "maf_with_groups_and_identifiers"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", help="comma-separated subset; default = all included")
    args = ap.parse_args()

    if args.accessions:
        accs = [a.strip() for a in args.accessions.split(",") if a.strip()]
    else:
        with open(REVIEW) as f:
            accs = [r["accession"] for r in csv.DictReader(f) if r["verdict"] == "include"]

    logger.info("Opening %d studies\n", len(accs))
    out_rows = []
    for acc in accs:
        try:
            files = list_files(acc)
        except Exception as exc:
            logger.error("%-12s file listing failed: %s", acc, exc)
            out_rows.append({"accession": acc, "verdict": "error",
                             "reason_code": f"file_listing_failed:{exc}"})
            continue

        names = [f["file"] for f in files if not f.get("directory")]
        mafs = sorted(n for n in names if n.startswith("m_") and n.endswith(".tsv"))
        sfiles = sorted(n for n in names if n.startswith("s_") and n.endswith(".txt"))
        assays = [n for n in names if n.startswith("a_") and n.endswith(".txt")]

        maf = None
        if mafs:
            # Largest MAF is the representative one; a study with several is
            # reported as such so the split is not silently collapsed.
            candidates = [inspect_maf(acc, n) for n in mafs]
            maf = max(candidates, key=lambda m: m["metabolites"])
            maf["metabolites_all"] = sum(c["metabolites"] for c in candidates)
            maf["identified_all"] = sum(c["identified"] for c in candidates)
            maf["with_db_id"] = sum(c["with_db_id"] for c in candidates)
            maf["with_name"] = sum(c["with_name"] for c in candidates)
            maf["sample_columns"] = sum(c["sample_columns"] for c in candidates)
            maf["declared_columns"] = sum(c["declared_columns"] for c in candidates)

        samples = inspect_samples(acc, sfiles[0]) if sfiles else {"n_samples": 0, "factors": {}}
        v, code = verdict(maf, samples)

        factor_desc = "; ".join(
            f"{k}={'/'.join(f'{lv}:{n}' for lv, n in sorted(levels.items()))}"
            for k, levels in samples["factors"].items()
        )
        row = {
            "accession": acc,
            "verdict": v,
            "reason_code": code,
            "n_maf": len(mafs),
            "n_assays": len(assays),
            "metabolites": maf["metabolites_all"] if maf else 0,
            "identified": maf["identified_all"] if maf else 0,
            "with_db_id": maf["with_db_id"] if maf else 0,
            "with_name": maf["with_name"] if maf else 0,
            "sample_columns": maf["sample_columns"] if maf else 0,
            "declared_columns": maf["declared_columns"] if maf else 0,
            "n_samples": samples["n_samples"],
            "factors": factor_desc,
        }
        out_rows.append(row)
        logger.info(
            "%-12s %-17s %5d metabolites (%5d named, %5d ChEBI) | %4d populated cols | "
            "%3d samples | %s",
            acc, v, row["metabolites"], row["with_name"], row["with_db_id"],
            row["sample_columns"], row["n_samples"],
            factor_desc[:62] or "(no factors)",
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "verdict", "reason_code", "n_maf", "n_assays",
              "metabolites", "identified", "with_db_id", "with_name",
              "sample_columns", "declared_columns", "n_samples", "factors"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    logger.info("\n-> %s", OUT.relative_to(REPO_ROOT))

    counts: dict[str, int] = {}
    for r in out_rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    logger.info("verdicts: %s", counts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
