"""Omics dataset ingestion pipeline step.

Queries public repositories for datasets matching the chronic pain search
criteria, filters by the configured thresholds, and writes a JSONL manifest
per modality to data/raw/{modality}/datasets.jsonl.

Available modalities: genomics, transcriptomics, proteomics, metabolomics, lipidomics
Run all five:  uv run python pipeline/03_ingest_omics.py --modality all
Run one:       uv run python pipeline/03_ingest_omics.py --modality transcriptomics
Dry run:       uv run python pipeline/03_ingest_omics.py --modality transcriptomics --dry-run
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

from cp_multiomics.ingest import IngesterFactory

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
DEFAULT_CONFIG = REPO_ROOT / "conf" / "ingest" / "default.yaml"
ALL_MODALITIES = ["genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"]


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def apply_filters(records: list, filters: dict) -> tuple[list, int]:
    """Drop records that fail min_samples or year_min thresholds."""
    min_samples = filters.get("min_samples", 0)
    year_min = int(filters.get("year_min", 0))
    kept, dropped = [], 0
    for r in records:
        if r.n_samples >= 0 and r.n_samples < min_samples:
            dropped += 1
            continue
        try:
            year = int(r.year[:4]) if r.year else 0
        except ValueError:
            year = 0
        if year_min and year and year < year_min:
            dropped += 1
            continue
        kept.append(r)
    return kept, dropped


def write_manifest(records: list, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for r in records:
            f.write(r.to_jsonl() + "\n")
    logger.info("Manifest: %d datasets → %s", len(records), out_path)


def run_modality(modality: str, cfg: dict, dry_run: bool) -> int:
    """Run ingestion for one modality; return number of datasets written."""
    search_cfg = cfg.get("search", {})
    # Mechanism phrases plus named conditions. The two are kept apart in the
    # config so that amending the condition list cannot move step 01's PRISMA
    # literature counts, and so the amendment is visible as its own change.
    terms = list(search_cfg.get("terms", ["chronic pain"]))
    terms += list(search_cfg.get("conditions", []))
    species = search_cfg.get("species", ["Homo sapiens"])
    filters = cfg.get("filters", {})
    output_cfg = cfg.get("output", {})
    hard_limit = output_cfg.get("hard_limit")          # None => paginate to exhaustion
    page_size = int(output_cfg.get("page_size", 100))
    raw_dir = REPO_ROOT / output_cfg.get("raw_dir", "data/raw")

    ingester = IngesterFactory(modality)

    if dry_run:
        query = ingester._build_query(terms, species)
        logger.info("[dry-run][%s] Would query: %s", modality, query[:200])
        return 0

    logger.info("=== Ingesting %s ===", modality)
    records = ingester.search_all(
        terms=terms,
        species_terms=species,
        filters=filters,
        hard_limit=int(hard_limit) if hard_limit else None,
        page_size=page_size,
    )

    records, dropped = apply_filters(records, filters)
    logger.info("[%s] After filters: %d kept, %d dropped", modality, len(records), dropped)

    out_path = raw_dir / modality / "datasets.jsonl"
    write_manifest(records, out_path)
    return len(records)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--modality",
        choices=ALL_MODALITIES + ["all"],
        required=True,
        help="Modality to ingest, or 'all' for all five",
    )
    p.add_argument(
        "--ingest-config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to ingest YAML config",
    )
    p.add_argument("--dry-run", action="store_true", help="Print queries without hitting APIs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.ingest_config)

    modalities = ALL_MODALITIES if args.modality == "all" else [args.modality]
    totals: dict[str, int] = {}
    failed: list[str] = []

    for modality in modalities:
        try:
            totals[modality] = run_modality(modality, cfg, dry_run=args.dry_run)
        except RuntimeError as exc:
            # Zero-recall is loud but must not abort the other modalities.
            logger.error("[%s] ingestion failed: %s", modality, exc)
            failed.append(modality)

    logger.info("=== Ingestion summary ===")
    for mod, count in totals.items():
        logger.info("  %-20s %d datasets", mod, count)
    for mod in failed:
        logger.error("  %-20s FAILED (zero recall or API error)", mod)

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
