"""GWAS Catalog FTP path resolution and harmonised-file download.

Harmonised sumstats live at a deterministic path keyed by accession, but the
filename varies, so the harmonised/ directory is listed and matched rather than
guessed. Some studies have no harmonised deposit (404) — that is a normal,
recorded exclusion, not an error.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from pathlib import Path

logger = logging.getLogger(__name__)


class FileTooLargeError(RuntimeError):
    """Raised when a harmonised file exceeds the configured size cap."""


def harmonised_dir(accession: str, ftp_base: str) -> str:
    """Return the harmonised/ directory URL for an accession."""
    n = int(accession.replace("GCST", ""))
    lo = ((n - 1) // 1000) * 1000 + 1
    hi = lo + 999
    return f"{ftp_base}/GCST{lo:06d}-GCST{hi:06d}/{accession}/harmonised/"


def find_harmonised_file(
    accession: str,
    ftp_base: str,
    list_dir: Callable[[str], list[str]],
) -> str | None:
    """Return the URL of the accession's *.h.tsv.gz, or None if absent.

    Prefers the canonical `{accession}.h.tsv.gz`; otherwise the first
    `*.h.tsv.gz` in the directory.
    """
    directory = harmonised_dir(accession, ftp_base)
    try:
        names = list_dir(directory)
    except Exception as exc:
        logger.warning("Listing failed for %s: %s", accession, exc)
        return None
    harmonised = [f for f in names if f.endswith(".h.tsv.gz")]
    if not harmonised:
        return None
    canonical = f"{accession}.h.tsv.gz"
    chosen = canonical if canonical in harmonised else sorted(harmonised)[0]
    return f"{directory}{chosen}"


def download_file(
    url: str,
    dest: Path,
    http_stream: Callable[[str], tuple[dict, Iterable[bytes]]],
    max_bytes: int,
) -> int:
    """Stream `url` to `dest`. Abort (no file written) if it exceeds max_bytes."""
    headers, chunks = http_stream(url)
    size = int(headers.get("Content-Length", 0) or 0)
    if size > max_bytes:
        raise FileTooLargeError(
            f"{url}: {size} bytes exceeds cap {max_bytes}"
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(dest, "wb") as f:
        for chunk in chunks:
            f.write(chunk)
            written += len(chunk)
    return written
