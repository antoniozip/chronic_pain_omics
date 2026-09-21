#!/usr/bin/env python
"""Download a published SPARC dataset from the Pennsieve Discover archive.

SPARC (sparc.science) hosts the NIH PRECISION Human Pain Network deposits,
which are CC-BY and need no access request. They are outside every repository
this project retrieves from -- GEO, PRIDE, MetaboLights and Metabolomics
Workbench -- so anything taken from SPARC is fetched here and screened in
`literature/prisma/sparc_candidates.csv`.

**Listing and downloading use two different hosts.** The Discover API browses a
version's file tree, but its `files/download` endpoint answers "The resource
requires authentication" even for a CC-BY dataset. The objects themselves are
public on the publish bucket, so the bytes come from there, keyed on the `s3://`
URI the listing already carries. Fetching a whole dataset through the API alone
fails with a 401 that reads like a permissions problem and is not one.

    python scripts/fetch_sparc_dataset.py 478 --dest data/raw/sparc/478

Sizes are checked against the listing, because a silent short read leaves a
truncated count matrix that parses.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import urllib.parse
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

API = "https://api.pennsieve.io/discover/datasets"
#: The publish bucket backing Discover. Public, and the only route that serves
#: a dataset's files without a token.
PUBLISH_BUCKET = "https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com"
TIMEOUT = "300"


def curl(url: str, out: Path | None = None) -> str:
    """GET a URL with curl, to a file when `out` is given."""
    cmd = ["curl", "-sfL", "--max-time", TIMEOUT, url]
    if out is not None:
        cmd += ["-o", str(out)]
    done = subprocess.run(cmd, capture_output=True, text=out is None)
    if done.returncode != 0:
        raise RuntimeError(f"curl failed ({done.returncode}) for {url}")
    return done.stdout if out is None else ""


def browse(dataset: int, version: int, path: str = "") -> list[dict]:
    """One directory level of a dataset version."""
    query = urllib.parse.urlencode({"path": path, "limit": 500})
    payload = json.loads(curl(f"{API}/{dataset}/versions/{version}/files/browse?{query}"))
    return payload.get("files") or []


def walk(dataset: int, version: int, path: str = ""):
    """Every file in the dataset, depth-first, as (relative path, size, uri)."""
    for entry in browse(dataset, version, path):
        child = f"{path}/{entry['name']}" if path else entry["name"]
        if entry.get("type") == "Directory":
            yield from walk(dataset, version, child)
        else:
            yield child, entry.get("size"), entry.get("uri", "")


def object_key(uri: str, dataset: int, relative: str) -> str:
    """The bucket key for a listed file.

    The listing's `s3://bucket/key` URI is authoritative: a key rebuilt from
    the dataset id and the relative path is right for the datasets seen so far
    and is a guess, so it is only the fallback.
    """
    if uri.startswith("s3://"):
        return uri.split("/", 3)[3]
    return f"{dataset}/{relative}"


def fetch(dataset: int, version: int, dest: Path) -> int:
    """Download the dataset, returning the number of files written."""
    written = 0
    for relative, size, uri in walk(dataset, version):
        out = dest / relative
        out.parent.mkdir(parents=True, exist_ok=True)
        curl(f"{PUBLISH_BUCKET}/{object_key(uri, dataset, relative)}", out)
        got = out.stat().st_size
        # A short read is the failure worth catching: a truncated CSV still
        # parses, with rows missing from the end and no error anywhere.
        if size is not None and got != size:
            logger.warning("%s: got %d bytes, listing says %d", relative, got, size)
        written += 1
        logger.info("%10d  %s", got, relative)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=int, help="SPARC dataset id, e.g. 478")
    parser.add_argument("--version", type=int, default=1)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args()

    n = fetch(args.dataset, args.version, args.dest)
    logger.info("wrote %d file(s) to %s", n, args.dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
