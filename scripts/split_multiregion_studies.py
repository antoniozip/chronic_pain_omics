#!/usr/bin/env python3
"""Split the two studies that deposit one count file per anatomical region.

GSE180627 and GSE197233 each deposit their counts region by region, and
``find_suppl_count_file()`` returns one, so each entered the pool contributing
a single region and discarding the rest -- the region chosen by filename sort
order. Merging the files instead is not an option: masking the region token
from the sample titles shows GSE180627's 96 samples are 24 animals x 4 regions
and GSE197233's 24 are 12 x 2, so a merged matrix would enter every animal
once per region.

This is the GSE241361 shape, and takes its treatment: a derived unit per
region, the parent retired to ``superseded_studies.csv``. Each derived unit
gets what ``run_transcriptomics_da()`` expects of an accession --

  {derived}_series_matrix.txt.gz   a copy of the parent's, so getGEO resolves
  {derived}/                       that region's count file, alone in it
  {derived}_groups.csv             the parent's curated arms, restricted

plus a record in ``data/interim/transcriptomics/harmonized.jsonl``, copied from
the parent, because ``05`` intersects ``--accessions`` with the manifest and a
derived id is in no manifest by construction. (GSE241361's splits are absent
from it, which is why that pair was produced outside the pipeline.)

-- so no pipeline code has to learn about derived ids.

Usage:
    python scripts/split_multiregion_studies.py [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_contrasts import CACHE, load_groups, parse_matrix, series_matrices

# parent -> (title separator, index of the region field, {region: [files]})
STUDIES = {
    "GSE180627": (
        "_", 2,
        {"PAG": ["GSE180627_Injury_PAG_TPM_counts.xlsx"],
         "PFC": ["GSE180627_Injury_PFC_TPM_counts.xlsx"],
         "S1":  ["GSE180627_Injury_S1_TPM_counts.xlsx"],
         "SC":  ["GSE180627_Injury_SC_TPM_counts.xlsx"]},
    ),
    # Sex splits the deposit but not the anatomy, and the two sexes are
    # different animals, so the region's two files merge into one unit.
    "GSE197233": (
        ".", 2,
        {"ACC":  ["GSE197233_Female-ACC-count.txt.gz",
                  "GSE197233_Male-ACC-count.txt.gz"],
         "mPFC": ["GSE197233_Female-mPFC-count.txt.gz",
                  "GSE197233_Male-mPFC-count.txt.gz"]},
    ),
}


def merge_tsv_gz(sources: list[Path], dest: Path) -> int:
    """Join count files on their first column; identical genes, disjoint samples."""
    merged: dict[str, list[str]] = defaultdict(list)
    header: list[str] = []
    for src in sources:
        with gzip.open(src, "rt") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        # GSE197233 deposits its column names with a trailing space
        # ("Female.Sham1 "), which no join to a sample title can survive.
        head, body = [c.strip() for c in rows[0]], rows[1:]
        header = header + head[1:] if header else list(head)
        for row in body:
            merged[row[0]].extend(row[1:])
    n_cols = len(header) - 1
    with gzip.open(dest, "wt") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for gene, vals in merged.items():
            if len(vals) == n_cols:          # drop genes absent from a file
                w.writerow([gene] + vals)
    return n_cols


MANIFEST = (Path(__file__).resolve().parents[1] / "data" / "interim"
            / "transcriptomics" / "harmonized.jsonl")


def register(parent: str, derived: str, region: str) -> bool:
    """Give the derived unit the parent's manifest record. Idempotent."""
    records = [json.loads(line) for line in MANIFEST.read_text().splitlines()
               if line.strip()]
    if any(r.get("accession") == derived for r in records):
        return False
    src = next((r for r in records if r.get("accession") == parent), None)
    if src is None:
        raise SystemExit(f"{parent} is not in {MANIFEST}")
    rec = dict(src)
    rec["accession"] = derived
    rec["derived_from"] = parent
    rec["region"] = region
    with open(MANIFEST, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for parent, (sep, idx, regions) in STUDIES.items():
        groups = load_groups(parent)
        matrix = series_matrices(parent)[0]
        gsms, _, titles = parse_matrix(matrix)
        by_region: dict[str, list[str]] = defaultdict(list)
        for gsm, title in zip(gsms, titles):
            parts = title.split(sep)
            if len(parts) > idx:
                by_region[parts[idx]].append(gsm)

        for region, files in regions.items():
            members = [g for g in by_region.get(region, []) if g in groups]
            arms = [(g, groups[g]) for g in members]
            n_case = sum(1 for _, a in arms if a == "case")
            derived = f"{parent}_{region}"
            print(f"{derived}: {n_case} case, {len(arms) - n_case} control "
                  f"from {len(files)} file(s)")
            if args.dry_run:
                continue
            if n_case < 3 or len(arms) - n_case < 3:
                print(f"  !! {derived}: an arm is below 3; skipped")
                continue

            out = CACHE / derived
            out.mkdir(exist_ok=True)
            srcs = [CACHE / parent / f for f in files]
            if len(srcs) == 1:
                shutil.copy(srcs[0], out / srcs[0].name)
            else:
                n = merge_tsv_gz(srcs, out / f"{derived}_raw_count_merged.txt.gz")
                print(f"  merged {len(srcs)} files -> {n} sample columns")
            shutil.copy(matrix, CACHE / f"{derived}_series_matrix.txt.gz")
            if register(parent, derived, region):
                print(f"  registered {derived} in harmonized.jsonl")
            with open(CACHE / f"{derived}_groups.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["sample_id", "group"])
                w.writeheader()
                w.writerows({"sample_id": g, "group": a} for g, a in arms)
    return 0


if __name__ == "__main__":
    sys.exit(main())
