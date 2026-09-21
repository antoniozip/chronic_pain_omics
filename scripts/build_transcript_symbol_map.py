#!/usr/bin/env python3
"""Build an Ensembl transcript -> gene symbol map for pipeline/05b.

Why this exists
---------------
``detect_id_space()`` routes human Ensembl *transcript* accessions to
``map_ensembl_transcript()``, which needs a lookup table. The two sources 05b
would normally reach cannot supply one here:

* biomaRt answers HTTP 404 from ``useMart()`` on this machine (see CLAUDE.md),
  and it is what the rat transcript map was originally built with.
* ``org.Hs.eg.db`` carries roughly 39,000 ``ENSEMBLTRANS`` keys against more
  than 250,000 human transcripts, and resolved 1 of 5 sampled ids from
  GSE153739/GSE153740. Below the match floor, so 05b correctly declines.

MyGene.info is already this project's ortholog service
(``src/cp_multiomics/ortholog/hcop.py``), it is reachable, and it indexes
``ensembl.transcript``. This script queries it once and writes the cache file
05b reads, in the same two-column form as ``rat_ensrnot2symbol.csv``.

The output is an input to a pipeline step, not a result: it lands in
``data/interim/`` and is re-usable by any later study on the same id space.

Usage
-----
    python scripts/build_transcript_symbol_map.py --species human
    python scripts/build_transcript_symbol_map.py --species human --dry-run
    python scripts/build_transcript_symbol_map.py --species human \
        --studies GSE153739 GSE153740
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
PER_STUDY = REPO_ROOT / "results" / "per_study" / "transcriptomics"
INTERIM = REPO_ROOT / "data" / "interim"

MYGENE_QUERY_URL = "https://mygene.info/v3/query"
BATCH_SIZE = 1000
RATE_DELAY = 0.34          # MyGene.info asks for <= 3 requests/second
MAX_RETRIES = 4

#: Transcript accession prefix and MyGene species tag, per species.
SPECIES: dict[str, dict[str, str]] = {
    "human": {"prefix": "ENST", "tag": "human", "cache": "human_enst2symbol.csv"},
    "mouse": {"prefix": "ENSMUST", "tag": "mouse", "cache": "mouse_ensmust2symbol.csv"},
    "rat": {"prefix": "ENSRNOT", "tag": "rat", "cache": "rat_ensrnot2symbol.csv"},
}

logger = logging.getLogger(__name__)


def collect_transcript_ids(prefix: str, studies: list[str] | None) -> list[str]:
    """Every transcript accession present in the per-study effect tables.

    Reads the *normalized* tables when they exist and the raw ones otherwise,
    which is the same preference ``06`` applies when it globs this directory.
    Version suffixes are stripped: Ensembl indexes accessions without them.
    """
    if not PER_STUDY.exists():
        raise SystemExit(f"no per-study directory at {PER_STUDY}")

    # One table per study, normalized preferred -- the same preference 06
    # applies when it globs this directory.
    chosen: dict[str, Path] = {}
    for path in sorted(PER_STUDY.glob("*_effects*.csv")):
        study = re.sub(r"_effects(_normalized)?\.csv$", "", path.name)
        if studies and study not in studies:
            continue
        if study not in chosen or path.name.endswith("_effects_normalized.csv"):
            chosen[study] = path

    found: set[str] = set()
    for study, path in sorted(chosen.items()):
        ids = pd.read_csv(path, usecols=["feature_id"], dtype=str)["feature_id"]
        hits = ids[ids.astype(str).str.startswith(prefix)]
        if len(hits):
            logger.info("  %s: %d %s ids", study, len(hits), prefix)
            found |= set(hits.astype(str).str.replace(r"\.\d+$", "", regex=True))
    return sorted(found)


def query_mygene(batch: list[str], tag: str) -> dict[str, str]:
    """One MyGene.info POST, retried with backoff. Returns only real hits."""
    payload = {
        "q": ",".join(batch),
        "scopes": "ensembl.transcript",
        "species": tag,
        "fields": "symbol",
        "size": "1",
    }
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.post(MYGENE_QUERY_URL, data=payload, timeout=120)
            resp.raise_for_status()
            hits = resp.json()
            break
        except Exception as exc:                      # noqa: BLE001 - retried
            wait = 2**attempt
            logger.warning("batch failed (%s); retrying in %ds", exc, wait)
            time.sleep(wait)
    else:
        logger.error("batch of %d abandoned after %d attempts",
                     len(batch), MAX_RETRIES)
        return {}

    out: dict[str, str] = {}
    for hit in hits:
        if hit.get("notfound"):
            continue
        symbol = hit.get("symbol")
        query = hit.get("query")
        # A transcript can hit more than one gene record; keep the first, which
        # is MyGene's best-scoring, and never overwrite it with a later tie.
        if symbol and query and query not in out:
            out[str(query)] = str(symbol)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", choices=sorted(SPECIES), default="human")
    parser.add_argument("--studies", nargs="*", default=None,
                        help="restrict to these accessions (default: all)")
    parser.add_argument("--dry-run", action="store_true",
                        help="report how many ids would be queried, then stop")
    parser.add_argument("--out-dir", type=Path, default=INTERIM)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    spec = SPECIES[args.species]
    cache_path = args.out_dir / spec["cache"]

    logger.info("collecting %s transcript ids from %s", spec["prefix"], PER_STUDY)
    ids = collect_transcript_ids(spec["prefix"], args.studies)
    if not ids:
        logger.info("no %s ids found; nothing to do", spec["prefix"])
        return 0
    logger.info("%d unique %s accessions", len(ids), spec["prefix"])

    # Resume rather than refetch: the file is the unit of work, not the run.
    known: dict[str, str] = {}
    if cache_path.exists():
        prev = pd.read_csv(cache_path, dtype=str)
        known = dict(zip(prev["ensembl_transcript_id"], prev["external_gene_name"]))
        logger.info("cache holds %d accessions already", len(known))
    todo = [i for i in ids if i not in known]
    logger.info("%d to query in %d batches", len(todo),
                (len(todo) + BATCH_SIZE - 1) // BATCH_SIZE)

    if args.dry_run:
        return 0

    resolved: dict[str, str] = {}
    for start in range(0, len(todo), BATCH_SIZE):
        batch = todo[start : start + BATCH_SIZE]
        resolved.update(query_mygene(batch, spec["tag"]))
        done = min(start + BATCH_SIZE, len(todo))
        logger.info("  %d/%d queried, %d resolved so far",
                    done, len(todo), len(resolved))
        if done < len(todo):
            time.sleep(RATE_DELAY)

    known.update(resolved)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame = (pd.DataFrame({"ensembl_transcript_id": list(known),
                           "external_gene_name": list(known.values())})
             .sort_values("ensembl_transcript_id"))
    frame.to_csv(cache_path, index=False, quoting=1)

    rate = 100.0 * sum(1 for i in ids if i in known) / len(ids)
    logger.info("wrote %d rows -> %s", len(frame), cache_path)
    logger.info("coverage of the ids seen in the corpus: %.1f%%", rate)
    return 0


if __name__ == "__main__":
    sys.exit(main())
