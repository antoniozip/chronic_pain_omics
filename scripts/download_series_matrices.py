#!/usr/bin/env python
"""Download GEO series matrix files for screened candidate accessions.

`auto_annotate_geo.py` derives case/control splits from
`data/raw/geo_cache/{accession}_series_matrix.txt.gz`, but nothing in the
repository fetched those files -- the original 110 were obtained outside the
pipeline, so re-running retrieval produced 165 candidates with no way to
annotate them. This closes that gap.

A series measured on several platforms publishes one matrix per platform
(`GSE12345-GPL570_series_matrix.txt.gz`). The directory listing is read rather
than a filename guessed, so multi-platform series are fetched completely
instead of silently yielding one platform's samples.

Existing files are skipped, so the script is resumable: GEO throttles, and a
run that dies partway must not restart from nothing.

Usage:
    python scripts/download_series_matrices.py --candidates <csv> [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
import time
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE = REPO_ROOT / "data" / "raw" / "geo_cache"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/matrix/"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)
SESSION = requests.Session()
TIMEOUT = 120


def stem_for(acc: str) -> str:
    """GSE67311 -> GSE67nnn, the FTP directory grouping GEO uses."""
    digits = acc[3:]
    return "GSE" + (digits[:-3] if len(digits) > 3 else "") + "nnn"


def matrix_names(acc: str) -> list[str]:
    url = BASE.format(stem=stem_for(acc), acc=acc)
    try:
        r = SESSION.get(url, timeout=TIMEOUT)
        r.raise_for_status()
    except Exception as exc:
        logger.warning("  %s: listing failed: %s", acc, exc)
        return []
    return sorted(set(re.findall(r'href="([^"]*_series_matrix\.txt\.gz)"', r.text)))


def fetch(acc: str, name: str) -> bool:
    out = CACHE / name
    if out.exists() and out.stat().st_size > 0:
        return True
    url = BASE.format(stem=stem_for(acc), acc=acc) + name
    try:
        with SESSION.get(url, timeout=TIMEOUT, stream=True) as r:
            r.raise_for_status()
            tmp = out.with_suffix(out.suffix + ".part")
            with tmp.open("wb") as fh:
                for chunk in r.iter_content(1 << 16):
                    fh.write(chunk)
            tmp.rename(out)          # rename last: a partial file must never
        return True                  # look like a complete one to the annotator
    except Exception as exc:
        logger.warning("  %s: %s failed: %s", acc, name, exc)
        return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    with args.candidates.open() as f:
        accs = [r["accession"] for r in csv.DictReader(f)
                if r["verdict"] == "pending"
                and r["retrieval"].startswith("amendment")]
    if args.limit:
        accs = accs[:args.limit]
    CACHE.mkdir(parents=True, exist_ok=True)

    todo = [a for a in accs if not list(CACHE.glob(f"{a}*_series_matrix.txt.gz"))]
    logger.info("%d candidates | %d already cached | %d to fetch",
                len(accs), len(accs) - len(todo), len(todo))

    ok = failed = 0
    for i, acc in enumerate(todo, 1):
        names = matrix_names(acc)
        if not names:
            failed += 1
            continue
        got = sum(fetch(acc, n) for n in names)
        ok += 1 if got else 0
        failed += 0 if got else 1
        if i % 20 == 0:
            logger.info("  %d/%d (%d ok, %d failed)", i, len(todo), ok, failed)
        time.sleep(0.34)             # GEO throttles sustained clients
    logger.info("\n%d fetched, %d failed", ok, failed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
