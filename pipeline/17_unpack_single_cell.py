"""Unpack the GEO _RAW.tar archives fetched by step 16.

Each archive holds one file per GSM, or a triplet of 10x files per GSM. They
are extracted into `data/raw/single_cell/{accession}/unpacked/` and a manifest
records which GSMs appeared, so the contrast audit has a list of samples that
actually carry data rather than one taken from the series metadata.

Only regular files are extracted, and any member whose path escapes the target
directory is refused: a tar from an external source should not be trusted to
stay inside the directory it is unpacked into.

    python pipeline/17_unpack_single_cell.py
    python pipeline/17_unpack_single_cell.py --accession GSE213216
"""

from __future__ import annotations

import argparse
import logging
import re
import tarfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_ROOT = REPO_ROOT / "data" / "raw" / "single_cell"
MANIFEST = RAW_ROOT / "unpack_manifest.csv"

GSM_RE = re.compile(r"(GSM\d+)")


def safe_members(tar: tarfile.TarFile, target: Path):
    """Yield the regular-file members that stay inside `target`."""
    for member in tar.getmembers():
        if not member.isfile():
            continue
        resolved = (target / member.name).resolve()
        if not resolved.is_relative_to(target.resolve()):
            logger.warning("  refusing member escaping the target: %s", member.name)
            continue
        yield member


def unpack(archive: Path, target: Path) -> list[str]:
    """Extract one archive. Returns the member names written."""
    target.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    with tarfile.open(archive) as tar:
        for member in safe_members(tar, target):
            dest = target / member.name
            # Re-extracting a 17 GB archive on every run is not free; a member
            # already present at the recorded size is left alone.
            if dest.exists() and dest.stat().st_size == member.size:
                written.append(member.name)
                continue
            tar.extract(member, path=target, filter="data")
            written.append(member.name)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", action="append")
    args = parser.parse_args()

    archives = sorted(RAW_ROOT.glob("*/[!.]*_RAW.tar"))
    if args.accession:
        archives = [a for a in archives if a.parent.name in args.accession]
    if not archives:
        logger.error("no _RAW.tar found under %s", RAW_ROOT)
        return 1

    rows: list[dict] = []
    for archive in archives:
        accession = archive.parent.name
        target = archive.parent / "unpacked"
        try:
            members = unpack(archive, target)
        except (tarfile.TarError, OSError) as exc:
            logger.error("%s: %s", accession, exc)
            rows.append({"accession": accession, "gsm": "", "n_files": 0,
                         "outcome": f"failed: {exc}"})
            continue

        gsms = sorted({m.group(1) for name in members
                       if (m := GSM_RE.search(name))})
        size = sum(f.stat().st_size for f in target.rglob("*") if f.is_file())
        logger.info("%s: %d member(s), %d GSM(s), %.2f GB",
                    accession, len(members), len(gsms), size / 1e9)
        for gsm in gsms or [""]:
            rows.append({"accession": accession, "gsm": gsm,
                         "n_files": sum(1 for n in members if gsm and gsm in n),
                         "outcome": "unpacked"})

    frame = pd.DataFrame(rows)
    frame.to_csv(MANIFEST, index=False)
    ok = frame[frame["outcome"] == "unpacked"]
    logger.info("%d archive(s), %d GSM(s) with files -> %s",
                frame["accession"].nunique(), len(ok[ok["gsm"] != ""]), MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
