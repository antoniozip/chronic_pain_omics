#!/usr/bin/env python3
"""Regenerate Table 4's two cross-species sensitivity arms.

Each arm asks what the rodent-to-human concordance becomes when the *human*
pool is restricted; the rodent pools are identical across all three arms, which
is what makes the comparison a statement about tissue compartment rather than
about pooling depth.

    published_pool   the human studies retrieved before the 2026-08-28 search
                     amendment, read from the screening sheet's `retrieval`
                     column rather than listed here
    tissue_matched   the human studies whose compartment is comparable to the
                     rodent pools -- nerve, spinal cord and blood -- read from
                     conf/analysis/human_tissue_compartments.csv

Until now no script did this. `07_cross_species.py` describes the arms running
06g into a work directory, and the sequence lived in someone's shell history,
so the arms sat three weeks stale behind a current whole-corpus arm twice: once
from 2026-08-29, and again after the 2026-09-12 corpus repair. A comparison
whose rows come from different corpora is worse than no comparison, because it
reads as a trend.

Every human study unit in the pool must carry a compartment. An unannotated one
is an error rather than a silent omission: it would otherwise fall out of
`tissue_matched` simply for being new, which is the failure this whole table
exists to detect.

Usage:
    python scripts/build_sensitivity_arms.py
    python scripts/build_sensitivity_arms.py --arm tissue_matched --dry-run
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
COMPARTMENTS = REPO_ROOT / "conf" / "analysis" / "human_tissue_compartments.csv"
CANDIDATES = REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv"


def human_units(meta_dir: Path, modality: str) -> list[str]:
    """Study units in the unrestricted human pool, as 06g wrote it."""
    pooled = meta_dir / modality / "by_species" / "Homo_sapiens_pooled.csv"
    if not pooled.exists():
        raise SystemExit(f"run 06g first: {pooled} is missing")
    d = pd.read_csv(pooled, usecols=["study_ids"])
    return sorted(d["study_ids"].astype(str).str.split(";")
                  .explode().str.strip().unique())


def _parent(unit: str) -> str:
    """A derived unit inherits its parent accession's curation.

    `GSE241361_DRG` is created by the analysis and appears in no curated sheet,
    so a lookup keyed on the accession answers "no match" for every split
    study -- which reads as "nothing to do" and is how five maps went quietly
    wrong on 2026-08-30.
    """
    return unit.split("_")[0]


def keep_published(units: list[str]) -> set[str]:
    d = pd.read_csv(CANDIDATES, dtype=str).set_index("accession")
    return {u for u in units
            if _parent(u) in d.index
            and d.loc[_parent(u), "retrieval"] == "original"}


def keep_tissue_matched(units: list[str]) -> set[str]:
    d = pd.read_csv(COMPARTMENTS, dtype=str).set_index("study_id")["compartment"]
    missing = [u for u in units if _parent(u) not in d.index]
    if missing:
        raise SystemExit(
            "these human study units carry no compartment in "
            f"{COMPARTMENTS.relative_to(REPO_ROOT)}: {', '.join(missing)}\n"
            "Annotate them; leaving them out would drop them from the "
            "tissue-matched arm for being new rather than for being a "
            "different compartment.")
    return {u for u in units if d.loc[_parent(u)] == "nerve_blood"}


ARMS = {"published_pool": keep_published, "tissue_matched": keep_tissue_matched}


def run_arm(name: str, modality: str, meta_dir: Path, out_root: Path,
            dry_run: bool) -> None:
    units = human_units(meta_dir, modality)
    keep = ARMS[name](units)
    if not keep:
        raise SystemExit(f"{name}: no human study units selected")
    excluded = sorted(set(units) - keep)
    print(f"[{name}] keeping {len(keep)} of {len(units)} human units: "
          f"{', '.join(sorted(keep))}")

    work = Path(tempfile.mkdtemp(prefix=f"sens_{name}_"))
    try:
        exclude_csv = work / "exclude.csv"
        exclude_csv.write_text("study_id\n" + "\n".join(excluded) + "\n")
        by_species = work / "meta" / modality / "by_species"
        by_species.mkdir(parents=True)
        out_dir = out_root / modality / "sensitivity" / name

        steps = [
            ["Rscript", "pipeline/06g_species_meta.R",
             "--modality", modality,
             "--exclude-studies", str(exclude_csv),
             "--out-dir", str(by_species)],
            [sys.executable, "pipeline/07_cross_species.py",
             "--modality", modality,
             "--meta-dir", str(work / "meta"),
             "--out-dir", str(out_dir)],
        ]
        for cmd in steps:
            print(f"[{name}] $ {' '.join(cmd)}")
            if dry_run:
                continue
            r = subprocess.run(cmd, cwd=REPO_ROOT)
            if r.returncode != 0:
                raise SystemExit(f"[{name}] failed: {' '.join(cmd)}")
        if not dry_run:
            summary = out_dir / "summary.csv"
            if not summary.exists():
                raise SystemExit(
                    f"[{name}] produced no {summary.name}. 07 exits 0 when it "
                    "finds nothing to compare, so an absent summary is how "
                    "these arms went stale while the chain reported success.")
            print(f"[{name}] wrote {summary.relative_to(REPO_ROOT)}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modality", default="transcriptomics")
    ap.add_argument("--arm", choices=sorted(ARMS), action="append",
                    help="repeatable; default is every arm")
    ap.add_argument("--meta-dir", type=Path,
                    default=REPO_ROOT / "results" / "meta")
    ap.add_argument("--out-root", type=Path,
                    default=REPO_ROOT / "results" / "cross_species")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for name in args.arm or sorted(ARMS):
        run_arm(name, args.modality, args.meta_dir, args.out_root, args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
