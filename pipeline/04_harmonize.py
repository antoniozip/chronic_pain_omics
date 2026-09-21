"""Harmonization pipeline step.

Reads raw ingest manifests (data/raw/{modality}/datasets.jsonl) produced by
step 03 and writes harmonized metadata tables to data/interim/{modality}/harmonized.jsonl.

Harmonization at this step is metadata-level:
  - Platform → controlled vocabulary (RNA-seq, GWAS, LC-MS, …)
  - Species → NCBI Taxonomy canonical names
  - ID space annotation (HGNC / rsID / UniProt / ChEBI / LIPIDMAPS)
  - Effect-size unit annotation (log2FC / beta / OR / SMD)
  - Ortholog flag for cross-species records
  - Quality flags (low_n, missing_year, …)

Full identifier remapping (gene symbol → ENSEMBL, rsID liftover, etc.) is
applied at the per-study DA stage (step 05) where raw count matrices are available.

Usage:
    uv run python pipeline/04_harmonize.py --modality transcriptomics
    uv run python pipeline/04_harmonize.py --modality all
    uv run python pipeline/04_harmonize.py --modality all --dry-run
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import yaml

from cp_multiomics.harmonize import HarmonizerFactory

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
DEFAULT_CONFIG = REPO_ROOT / "conf" / "ingest" / "default.yaml"
ALL_MODALITIES = ["genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"]


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_raw_records(raw_dir: Path, modality: str) -> list[dict]:
    jsonl_path = raw_dir / modality / "datasets.jsonl"
    if not jsonl_path.exists():
        logger.warning("No raw manifest for %s at %s — skipping", modality, jsonl_path)
        return []
    records: list[dict] = []
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    logger.info("[%s] Loaded %d raw records from %s", modality, len(records), jsonl_path)
    return records


def run_modality(modality: str, cfg: dict, raw_dir: Path, interim_dir: Path, dry_run: bool) -> int:
    raw_records = load_raw_records(raw_dir, modality)
    if not raw_records:
        return 0

    if dry_run:
        logger.info("[dry-run][%s] Would harmonize %d records", modality, len(raw_records))
        return 0

    harmonizer = HarmonizerFactory(modality)
    harmonized = harmonizer.harmonize(raw_records, cfg)

    out_dir = interim_dir / modality
    harmonizer.save(harmonized, out_dir)

    # Log quality summary
    flagged = [r for r in harmonized if r.quality_flags]
    ortholog_needed = sum(1 for r in harmonized if r.ortholog_needed)
    logger.info(
        "[%s] %d harmonized | %d quality-flagged | %d need ortholog mapping",
        modality, len(harmonized), len(flagged), ortholog_needed,
    )
    return len(harmonized)


def print_summary(modality: str, interim_dir: Path) -> None:
    jsonl_path = interim_dir / modality / "harmonized.jsonl"
    if not jsonl_path.exists():
        return
    records = [json.loads(line) for line in jsonl_path.read_text().splitlines() if line.strip()]

    platforms = {}
    id_spaces = {}
    for r in records:
        p = r.get("platform_canonical", "unknown")
        platforms[p] = platforms.get(p, 0) + 1
        i = r.get("id_space", "unknown")
        id_spaces[i] = id_spaces.get(i, 0) + 1

    print(f"\n  [{modality}] {len(records)} datasets")
    for p, n in sorted(platforms.items(), key=lambda x: -x[1]):
        print(f"    platform  {p:<30} {n:>4}")
    for i, n in sorted(id_spaces.items(), key=lambda x: -x[1]):
        print(f"    id_space  {i:<30} {n:>4}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--modality",
        choices=ALL_MODALITIES + ["all"],
        required=True,
    )
    p.add_argument("--ingest-config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--raw-dir", type=Path, default=REPO_ROOT / "data" / "raw")
    p.add_argument("--interim-dir", type=Path, default=REPO_ROOT / "data" / "interim")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.ingest_config)
    modalities = ALL_MODALITIES if args.modality == "all" else [args.modality]

    totals: dict[str, int] = {}
    for modality in modalities:
        totals[modality] = run_modality(
            modality, cfg, args.raw_dir, args.interim_dir, args.dry_run
        )

    logger.info("=== Harmonization summary ===")
    if not args.dry_run:
        print("\n=== Harmonized dataset breakdown ===")
        for modality in modalities:
            print_summary(modality, args.interim_dir)
    for mod, n in totals.items():
        logger.info("  %-20s %d records", mod, n)


if __name__ == "__main__":
    main()
