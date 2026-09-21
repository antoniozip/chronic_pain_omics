#!/usr/bin/env python
"""Build the smallest results archive that satisfies every manuscript claim.

The CI claims job needs the pipeline outputs, but `results/` is ~12 GB, well
past the 10 GB GitHub Actions allows a repository in cache. It does not need
all of it: as of 2026-08-28 the claim registry reads 40 files totalling 140 MB,
which compress to 45 MB.

The file list is derived by running the registry and recording what it opens,
so it follows the claims automatically. Adding a claim that reads a new results
file changes this archive without anyone maintaining a list — which also means
the counts above drift (16 files, then 23, now 40). This script prints the
current ones; prefer them to any figure written down elsewhere.

Usage:
    python scripts/build_claims_archive.py                    # -> results-claims.tar.gz
    python scripts/build_claims_archive.py --out /tmp/r.tgz
"""

from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tests"))

import test_manuscript_claims as tmc  # noqa: E402


def required_files() -> set[Path]:
    """Every results file the claim registry opens."""
    seen: set[Path] = set()
    original = tmc._read

    def spy(path, **kwargs):
        seen.add(Path(path))
        return original(path, **kwargs)

    tmc._read = spy
    try:
        for claim in [*tmc.CLAIMS, *tmc.GENE_CLAIMS]:
            try:
                claim.expected()
            except Exception:
                # A claim whose input is missing still records the path it
                # wanted, which is the point of this pass.
                pass
    finally:
        tmc._read = original
    return seen


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "results-claims.tar.gz")
    args = parser.parse_args()

    needed = sorted(required_files())
    missing = [p for p in needed if not p.exists()]
    if missing:
        print("cannot build: these inputs have not been generated:")
        for p in missing:
            print(f"  - {p.relative_to(REPO_ROOT)}")
        print("\nRun the pipeline first; a partial archive would let CI skip claims.")
        return 1

    total = 0
    with tarfile.open(args.out, "w:gz") as tar:
        for path in needed:
            total += path.stat().st_size
            tar.add(path, arcname=str(path.relative_to(REPO_ROOT)))

    size = args.out.stat().st_size
    print(f"{len(needed)} files, {total / 1e6:.1f} MB -> {args.out} ({size / 1e6:.1f} MB)")
    print(
        "\nUpload it somewhere the runner can reach over HTTPS, then run the CI "
        "workflow manually with results_url set to that address."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
