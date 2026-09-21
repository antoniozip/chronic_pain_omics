"""Capture external-tool versions and per-run provenance.

METAL and MAGMA are pinned binaries, so neither `uv.lock` nor `renv.lock`
records what actually produced a result. Without that, a genomics result cannot
be reproduced from the repository alone. Every stage that shells out to one of
them writes a `run_provenance.json` beside its outputs.

A version probe must never break an analysis, so every failure mode here
degrades to "unknown" rather than raising.
"""

from __future__ import annotations

import json
import logging
import subprocess
from datetime import UTC, datetime
from pathlib import Path

logger = logging.getLogger(__name__)

UNKNOWN = "unknown"
_PROBE_TIMEOUT_S = 15.0

# Errors that mean "this binary did not tell us its version", never a reason to
# abort the surrounding analysis.
_PROBE_ERRORS = (
    OSError,                       # missing, not executable, not a binary
    subprocess.TimeoutExpired,     # waited on stdin, or hung
    ValueError,
)


def probe_version(
    binary: Path,
    args: list[str] | None = None,
    stdin: str = "",
    timeout: float = _PROBE_TIMEOUT_S,
) -> str:
    """Return the first banner line the binary prints, or "unknown".

    `stdin` is fed as empty by default so a tool that reads a command script
    from standard input (METAL) sees EOF and exits instead of blocking.
    """
    try:
        result = subprocess.run(
            [str(binary), *(args or [])],
            input=stdin, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except _PROBE_ERRORS as exc:
        logger.warning("could not probe %s version: %s", binary, exc)
        return UNKNOWN

    if result.returncode != 0 and not result.stdout.strip():
        logger.warning("%s exited %d without a version banner",
                       binary, result.returncode)
        return UNKNOWN

    for line in (result.stdout + result.stderr).splitlines():
        text = line.strip()
        if text:
            return text
    return UNKNOWN


def git_revision(repo_root: Path) -> str:
    """Return the current commit SHA, or "unknown" outside a git checkout."""
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=_PROBE_TIMEOUT_S, check=False,
        )
    except _PROBE_ERRORS as exc:
        logger.warning("could not read git revision: %s", exc)
        return UNKNOWN
    return result.stdout.strip() or UNKNOWN


def write_run_provenance(
    out_dir: Path,
    tool: str,
    binary: Path,
    version: str,
    config_path: Path | None = None,
    repo_root: Path | None = None,
    extra: dict | None = None,
) -> Path:
    """Write run_provenance.json into `out_dir`; return its path."""
    record = {
        "tool": tool,
        "binary": str(binary),
        "version": version,
        "config": str(config_path) if config_path else None,
        "git_revision": git_revision(repo_root or Path.cwd()),
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    if extra:
        record.update(extra)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "run_provenance.json"
    out_path.write_text(json.dumps(record, indent=2) + "\n")
    logger.info("%s provenance (%s) -> %s", tool, version, out_path)
    return out_path
