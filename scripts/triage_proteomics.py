#!/usr/bin/env python
"""Phase 2 triage: which discovered proteomics datasets are actually usable?

Discovery (03) says a dataset is about pain. It says nothing about whether the
submitters deposited anything we can compute an effect size from. PXD013362 was
usable because it shipped a processed quantification table with group labels;
a deposit of .raw files alone is a paper we can read, not data we can pool.

Classifies every harmonized record by what PRIDE holds for it:

  processed_table  a protein/peptide quantification matrix is present
  raw_only         only instrument or peak files
  ineligible       wrong organism or not a pain phenotype
  unresolved       PRIDE file listing failed

`processed_table` is a *candidate*, not a verdict: only opening the file shows
whether it carries per-sample values and recoverable group labels. The column
is named accordingly.

Writes literature/prisma/proteomics_triage.csv with a reason code per record,
per the PRISMA discipline in CLAUDE.md.
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "data" / "interim" / "proteomics" / "harmonized.jsonl"
OUT = REPO_ROOT / "literature" / "prisma" / "proteomics_triage.csv"
FILES_URL = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{}/files"

ELIGIBLE_SPECIES = {"Homo sapiens", "Mus musculus", "Rattus norvegicus"}

# Names that a protein- or peptide-level quantification table tends to carry.
# MaxQuant (proteinGroups/peptides), FragPipe (combined_protein), DIA-NN
# (report), plus the generic spreadsheet a submitter attaches as supplementary.
QUANT_PATTERNS = re.compile(
    r"(proteingroups|combined_protein|combined_peptide|report\.(tsv|pr_matrix)"
    r"|abundance|quantif|quantit|_quant|matrix|intensit|expression|lfq"
    r"|supplementary|supplemental|table)",
    re.I,
)
TABLE_EXT = (".txt", ".tsv", ".csv", ".xlsx", ".xls")
RAW_CATEGORIES = {"RAW", "PEAK"}


def file_listing(accession: str) -> list[dict] | None:
    for attempt in range(2):
        try:
            r = requests.get(FILES_URL.format(accession), timeout=60)
            r.raise_for_status()
            data = r.json()
            if isinstance(data, dict):
                data = data.get("_embedded", {}).get("files", data.get("files", []))
            return data if isinstance(data, list) else []
        except Exception:
            if attempt == 0:
                time.sleep(2)
                continue
            return None
    return None


def classify(record: dict, files: list[dict] | None) -> tuple[str, str]:
    species = set(record.get("species_canonical") or [])
    if not species & ELIGIBLE_SPECIES:
        return "ineligible", f"organism_not_eligible:{'|'.join(sorted(species)) or 'none'}"

    if files is None:
        return "unresolved", "pride_file_listing_failed"
    if not files:
        return "unresolved", "pride_reported_no_files"

    def category(f: dict) -> str:
        cat = f.get("fileCategory")
        return str(cat.get("value", "")) if isinstance(cat, dict) else str(cat or "")

    names = [(f.get("fileName") or "") for f in files]
    cats = {category(f) for f in files}

    hits = [n for n in names
            if n.lower().endswith(TABLE_EXT) and QUANT_PATTERNS.search(n)]
    if hits:
        return "processed_table", f"candidate_table:{hits[0][:60]}"

    if cats and all(c.upper() in RAW_CATEGORIES or not c for c in cats):
        return "raw_only", "instrument_or_peak_files_only"
    return "raw_only", "no_recognised_quantification_table"


def main() -> int:
    if not MANIFEST.exists():
        print(f"no manifest at {MANIFEST}; run 03 and 04 first")
        return 1
    records = [json.loads(line) for line in MANIFEST.read_text().splitlines() if line.strip()]
    print(f"triaging {len(records)} datasets against PRIDE ...")

    rows = []
    for i, rec in enumerate(records, 1):
        acc = rec["accession"]
        files = file_listing(acc)
        verdict, reason = classify(rec, files)
        rows.append({
            "accession": acc,
            "verdict": verdict,
            "reason_code": reason,
            "species": ";".join(sorted(rec.get("species_canonical") or [])),
            "has_human": rec.get("has_human"),
            "has_animal": rec.get("has_animal"),
            "n_files": 0 if files is None else len(files),
            "title": (rec.get("title") or "")[:120],
        })
        print(f"  [{i:2d}/{len(records)}] {acc}  {verdict:16s} {reason[:52]}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"\n-> {OUT.relative_to(REPO_ROOT)}")
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    for k in ("processed_table", "raw_only", "ineligible", "unresolved"):
        if k in counts:
            print(f"  {k:16s} {counts[k]}")

    usable = [r for r in rows if r["verdict"] == "processed_table"]
    print(f"\ncandidates: {len(usable)} "
          f"({sum(1 for r in usable if r['has_human'])} human, "
          f"{sum(1 for r in usable if r['has_animal'])} animal)")
    print("meta.min_studies is 3, so pooling needs at least 3 that survive "
          "manual inspection of the table itself.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
