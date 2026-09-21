#!/usr/bin/env python
"""Render the cover letter to .docx for submission.

The journal receives a Word file, not Markdown. This is the fourth Word file
the editor opens, and it goes through the same reference document as the
manuscript and the supplementary -- pandoc's stock styling
would otherwise leave it the only one not in Times.

The Markdown stays authoritative: it is what the claim-free prose is edited in,
and this script only renders it. Regenerate rather than edit the .docx.

Build-time only: pypandoc is not a project dependency, so run this with the
system python3 rather than .venv/bin/python.

Also renders the PRISMA 2020 checklist, which is Markdown for the same reason.

Usage:
    python3 scripts/build_cover_letter_docx.py
    python3 scripts/build_cover_letter_docx.py \\
        --source manuscript/prisma_2020_checklist.md \\
        --out manuscript/prisma_2020_checklist.docx
"""

from __future__ import annotations

import argparse
import logging
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import docx_reference  # noqa: E402
import pypandoc  # noqa: E402

SOURCE = REPO_ROOT / "manuscript" / "cover_letter.md"
OUT = REPO_ROOT / "manuscript" / "cover_letter.docx"

logger = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if not args.source.exists():
        raise SystemExit(f"missing source: {args.source}")

    text = args.source.read_text()
    with tempfile.TemporaryDirectory() as tmp:
        reference = docx_reference.build_reference(Path(tmp) / "reference.docx")
        pypandoc.convert_text(
            text, "docx", format="markdown",
            outputfile=str(args.out),
            extra_args=["--reference-doc", str(reference)],
        )
    logger.info("wrote %s (%.1f kB)",
                args.out.relative_to(REPO_ROOT), args.out.stat().st_size / 1e3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
