"""Stage 1 genomics: resolve, dedup, fetch, and harmonize GWAS summary statistics.

Retrieve + harmonize ONLY — no METAL/MAGMA (deferred). Applies the §7.7
cohort-independence rule BEFORE download, so only poolable studies are fetched.
Every study's disposition is recorded in a manifest (PRISMA discipline).

Usage:
    uv run python pipeline/04b_gwas_sumstats.py --dry-run   # resolve+dedup+manifest, no download
    uv run python pipeline/04b_gwas_sumstats.py --limit 2   # smoke test: fetch first 2 kept studies
    uv run python pipeline/04b_gwas_sumstats.py             # full fetch (~13 GB)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import requests
import yaml

from cp_multiomics.genomics.dedup import apply_presumed_cohort, dedup_studies
from cp_multiomics.genomics.ftp import (
    FileTooLargeError,
    download_file,
    find_harmonised_file,
)
from cp_multiomics.genomics.harmonize import harmonize_sumstats
from cp_multiomics.genomics.resolve import resolve_sumstats_studies

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
DEFAULT_CONFIG = REPO_ROOT / "conf" / "genomics" / "sumstats.yaml"

MANIFEST_COLUMNS = [
    "accession", "trait", "phenotype", "cohort", "n", "is_burden",
    "disposition", "reason", "harmonised_url", "local_path", "n_variants",
]


def _default_http_get(url: str, params: dict | None = None) -> dict[str, Any]:
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    return resp.json()


def _default_list_dir(url: str) -> list[str]:
    import re
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        return []
    return [
        f for f in re.findall(r'href="([^"]+)"', resp.text)
        if not f.startswith("/") and not f.startswith("?") and f != "../"
    ]


def _default_http_stream(url: str) -> tuple[dict, Iterable[bytes]]:
    resp = requests.get(url, stream=True, timeout=300)
    resp.raise_for_status()
    return dict(resp.headers), resp.iter_content(chunk_size=1 << 20)


def _verification_path(raw_path: Path) -> Path:
    """Sidecar recording that a cached raw file is a complete download."""
    return raw_path.with_suffix(raw_path.suffix + ".ok")


def write_verification(raw_path: Path) -> Path:
    """Write `{bytes, sha256}` for a freshly downloaded raw file."""
    digest = hashlib.sha256()
    with open(raw_path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    payload = {"bytes": raw_path.stat().st_size, "sha256": digest.hexdigest()}
    marker = _verification_path(raw_path)
    tmp = marker.with_suffix(marker.suffix + ".tmp")
    tmp.write_text(json.dumps(payload), encoding="utf-8")
    os.replace(tmp, marker)
    return marker


def raw_is_verified(raw_path: Path) -> bool:
    """True only when a prior download recorded matching verification metadata.

    A non-empty file is not proof of a complete download: an interrupted stream
    leaves a short but syntactically valid GWAS-SSF file, and re-runs used to
    harmonize it as if it were the whole deposit.
    """
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        return False
    marker = _verification_path(raw_path)
    if not marker.exists():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("bytes") == raw_path.stat().st_size


def build_manifest(kept, dropped, fetch_results: dict[str, dict]) -> list[dict]:
    """Assemble manifest rows. fetch_results maps accession -> partial row overrides."""
    rows: list[dict] = []
    for study, reason in dropped:
        rows.append({
            "accession": study.accession, "trait": study.trait,
            "phenotype": study.phenotype, "cohort": study.cohort, "n": study.n,
            "is_burden": study.is_burden, "disposition": "excluded",
            "reason": reason, "harmonised_url": "", "local_path": "", "n_variants": "",
        })
    for study in kept:
        base = {
            "accession": study.accession, "trait": study.trait,
            "phenotype": study.phenotype, "cohort": study.cohort, "n": study.n,
            "is_burden": study.is_burden, "disposition": "excluded",
            "reason": "", "harmonised_url": "", "local_path": "", "n_variants": "",
        }
        base.update(fetch_results.get(study.accession, {}))
        rows.append(base)
    return rows


def run(config, http_get=None, list_dir=None, http_stream=None,
        dry_run=False, limit=None) -> Path:
    http_get = http_get or _default_http_get
    list_dir = list_dir or _default_list_dir
    http_stream = http_stream or _default_http_stream

    out = config["output"]
    raw_dir = Path(out["raw_dir"])
    interim_dir = Path(out["interim_dir"])
    manifest_path = Path(out["manifest"])
    max_bytes = int(out.get("max_file_mb", 1500)) * 1_000_000

    studies = resolve_sumstats_studies(config["efo_traits"], http_get)
    studies = apply_presumed_cohort(studies)
    kept, dropped = dedup_studies(studies)
    logger.info("Resolved %d, kept %d after §7.7 dedup", len(studies), len(kept))

    to_fetch = kept if limit is None else kept[:limit]
    fetch_results: dict[str, dict] = {}

    for study in to_fetch:
        if dry_run:
            fetch_results[study.accession] = {"reason": "dry_run"}
            continue
        url = find_harmonised_file(study.accession, config["ftp_base"], list_dir)
        if url is None:
            fetch_results[study.accession] = {"reason": "no_harmonised_file"}
            continue
        raw_path = raw_dir / f"{study.accession}.h.tsv.gz"
        if raw_is_verified(raw_path):
            logger.info("[%s] using verified cached raw, skipping download",
                        study.accession)
        else:
            if raw_path.exists():
                logger.warning(
                    "[%s] cached raw is truncated or has no verification marker — "
                    "re-downloading", study.accession)
            try:
                download_file(url, raw_path, http_stream, max_bytes)
            except (FileTooLargeError, requests.RequestException, OSError) as exc:
                logger.error("[%s] download failed: %s", study.accession, exc)
                fetch_results[study.accession] = {
                    "reason": "download_failed", "harmonised_url": url}
                continue
            write_verification(raw_path)
        interim_path = interim_dir / f"{study.accession}.tsv"
        try:
            # A malformed real-world GWAS-SSF file (unexpected column set) is a
            # documented, expected failure mode (harmonize.py raises ValueError).
            # One bad file anywhere in a ~40-study run must never abort the batch
            # and destroy provenance for every other study, so this is caught
            # narrowly around just the harmonize call and recorded per-study.
            n_variants = harmonize_sumstats(raw_path, study, interim_path)
        except Exception as exc:
            logger.error("[%s] harmonize failed: %s", study.accession, exc)
            fetch_results[study.accession] = {
                "reason": "harmonize_failed", "harmonised_url": url}
            continue
        fetch_results[study.accession] = {
            "disposition": "fetched", "reason": "", "harmonised_url": url,
            "local_path": str(interim_path), "n_variants": str(n_variants),
        }

    if limit is not None:
        for study in kept[limit:]:
            fetch_results.setdefault(study.accession, {"reason": "not_attempted"})

    rows = build_manifest(kept, dropped, fetch_results)
    # A dry run is a cheap resolve+dedup probe, not a fetch: it must not
    # overwrite the PRISMA manifest with every study marked "excluded".
    write_path = (manifest_path.with_suffix(".dryrun.csv")
                  if dry_run else manifest_path)
    write_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = write_path.with_suffix(write_path.suffix + ".tmp")
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp_path, write_path)
    logger.info("Manifest written -> %s (%d rows)", write_path, len(rows))
    return write_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--dry-run", action="store_true",
                   help="Resolve + dedup + write manifest without downloading")
    p.add_argument("--limit", type=int, default=None,
                   help="Fetch only the first N kept studies (smoke test)")
    args = p.parse_args()
    with open(args.config) as f:
        config = yaml.safe_load(f)
    run(config, dry_run=args.dry_run, limit=args.limit)


if __name__ == "__main__":
    main()
