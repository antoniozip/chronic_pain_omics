"""Retrieve the single-cell and single-nucleus series into data/raw/single_cell/.

Fetches every supplementary file GEO publishes for the accessions in
`conf/analysis/single_cell_studies.csv`, one series to a directory, and writes
a manifest of what arrived.

Resumable by design. These are 10x matrix archives and a run over eighteen
series is long enough that it will be interrupted: a file whose local size
already equals the remote Content-Length is skipped, and a partial file is
re-fetched whole rather than appended to, since a truncated .tar that looks
complete is worse than one that is obviously missing.

    python pipeline/16_download_single_cell.py                 # all 18
    python pipeline/16_download_single_cell.py --accession GSE253345
    python pipeline/16_download_single_cell.py --dry-run       # sizes only
"""

from __future__ import annotations

import argparse
import logging
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
REGISTRY = REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv"
OUT_ROOT = REPO_ROOT / "data" / "raw" / "single_cell"
MANIFEST = OUT_ROOT / "download_manifest.csv"

USER_AGENT = "CP-multiomics/1.0 (academic meta-analysis)"
#: NCBI asks for restraint on its FTP endpoints; one series at a time with a
#: pause between requests is well inside what they publish as acceptable.
POLITE_DELAY_S = 0.5
CHUNK = 1 << 20


def series_url(accession: str) -> str:
    """GEO puts each series under a directory named for its thousand-block."""
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/"
            f"{accession}/suppl/")


def list_files(accession: str) -> list[str]:
    req = urllib.request.Request(series_url(accession), headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as response:
        html = response.read().decode(errors="replace")
    # The listing links every file relative to the directory; absolute hrefs
    # are the parent-directory and sort links, which are not files.
    return sorted({f for f in re.findall(r'href="([^"?/][^"]*)"', html)
                   if not f.startswith("http")})


def remote_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            length = response.headers.get("Content-Length")
        return int(length) if length else None
    except (urllib.error.URLError, ValueError):
        return None


def fetch(url: str, dest: Path, expected: int | None) -> tuple[str, int]:
    """Download one file. Returns (outcome, bytes on disk)."""
    if dest.exists() and expected is not None and dest.stat().st_size == expected:
        return "skipped", dest.stat().st_size
    if dest.exists():
        # Partial or size-unknown: start again rather than trust the remains.
        dest.unlink()

    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as response, tmp.open("wb") as out:
        while chunk := response.read(CHUNK):
            out.write(chunk)

    size = tmp.stat().st_size
    if expected is not None and size != expected:
        tmp.unlink()
        raise OSError(f"{dest.name}: got {size} bytes, expected {expected}")
    tmp.rename(dest)
    return "downloaded", size


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", action="append",
                        help="restrict to these accessions (repeatable)")
    parser.add_argument("--dry-run", action="store_true",
                        help="list files and total size without downloading")
    args = parser.parse_args()

    registry = pd.read_csv(REGISTRY, dtype=str)
    if args.accession:
        registry = registry[registry["accession"].isin(args.accession)]
        if registry.empty:
            logger.error("no registry entry matches %s", args.accession)
            return 2

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    grand_total = 0

    for _, study in registry.iterrows():
        accession = study["accession"]
        try:
            files = list_files(accession)
        except Exception as exc:                       # noqa: BLE001 - reported, not raised
            logger.error("%s: cannot list supplementary files (%s)", accession, exc)
            rows.append({"accession": accession, "file": "", "bytes": 0,
                         "outcome": f"listing_failed: {exc}"})
            continue

        target = OUT_ROOT / accession
        target.mkdir(parents=True, exist_ok=True)
        series_total = 0

        for name in files:
            url = series_url(accession) + name
            expected = remote_size(url)
            if args.dry_run:
                series_total += expected or 0
                rows.append({"accession": accession, "file": name,
                             "bytes": expected or 0, "outcome": "dry_run"})
                time.sleep(POLITE_DELAY_S)
                continue
            try:
                outcome, size = fetch(url, target / name, expected)
            except Exception as exc:                   # noqa: BLE001
                logger.error("  %s/%s failed: %s", accession, name, exc)
                rows.append({"accession": accession, "file": name, "bytes": 0,
                             "outcome": f"failed: {exc}"})
                continue
            series_total += size
            rows.append({"accession": accession, "file": name, "bytes": size,
                         "outcome": outcome})
            logger.info("  %s/%s  %.1f MB  %s", accession, name, size / 1e6, outcome)
            time.sleep(POLITE_DELAY_S)

        grand_total += series_total
        logger.info("%s: %d file(s), %.2f GB", accession, len(files), series_total / 1e9)

    pd.DataFrame(rows).to_csv(MANIFEST, index=False)
    logger.info("%d file(s) over %d series, %.2f GB total -> %s",
                len(rows), registry["accession"].nunique(), grand_total / 1e9, MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
