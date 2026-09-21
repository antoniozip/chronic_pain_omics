"""Rebuild effects_normalized_all.csv from the per-study normalized tables.

Step 05b performs two jobs: it maps platform identifiers to gene symbols, and
it concatenates the results. Re-running it to change only the second job is
not currently safe - the Bioconductor annotation packages it needs
(rat2302.db, rgu34a.db) are absent from renv.lock, so a full re-run silently
degrades nine studies to their native probe ids instead of gene symbols.

This script instead filters the existing aggregate in place, dropping only
the rows of superseded studies listed in conf/analysis/superseded_studies.csv.
Every retained row is preserved byte-for-byte. It deliberately does not
re-derive normalization, and it is not a substitute for re-running 05b once
the missing annotation packages are pinned.

    uv run python scripts/rebuild_normalized_aggregate.py
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DIR = REPO_ROOT / "results" / "per_study" / "transcriptomics"
SUPERSEDED_CSV = REPO_ROOT / "conf" / "analysis" / "superseded_studies.csv"


def load_superseded(path: Path) -> set[str]:
    """Study ids excluded because a finer-grained re-analysis supersedes them."""
    if not path.exists():
        logger.warning("no superseded registry at %s", path)
        return set()
    df = pd.read_csv(path)
    return set(df["study_id"].astype(str)) if "study_id" in df.columns else set()


def rebuild(per_study_dir: Path, superseded_csv: Path) -> Path:
    superseded = load_superseded(superseded_csv)
    out_path = per_study_dir / "effects_normalized_all.csv"
    if not out_path.exists():
        raise FileNotFoundError(f"no aggregate to filter at {out_path}")

    combined = pd.read_csv(out_path, low_memory=False)
    before_rows, before_studies = len(combined), combined["study_id"].nunique()

    present = superseded & set(combined["study_id"].astype(str))
    for study_id in sorted(present):
        n = int((combined["study_id"] == study_id).sum())
        logger.info("excluding %s: superseded, %d rows (see %s)",
                    study_id, n, superseded_csv.name)
    for study_id in sorted(superseded - present):
        logger.info("%s listed as superseded but not present", study_id)

    combined = combined[~combined["study_id"].astype(str).isin(superseded)]
    combined.to_csv(out_path, index=False)
    logger.info("%d -> %d rows, %d -> %d studies", before_rows, len(combined),
                before_studies, combined["study_id"].nunique())
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-study-dir", type=Path, default=DEFAULT_DIR)
    p.add_argument("--superseded-csv", type=Path, default=SUPERSEDED_CSV)
    args = p.parse_args()
    rebuild(args.per_study_dir, args.superseded_csv)


if __name__ == "__main__":
    main()
