"""Render every manuscript figure to SVG alongside its PDF.

Journals and co-authors routinely ask for vector sources that can be edited
without a LaTeX toolchain. PyMuPDF converts each PDF page to true vector SVG,
so nothing is rasterized.

Multi-page figures (the per-feature forest plots) produce one SVG per page,
named `<stem>_p01.svg`; single-page figures keep the plain stem.

By default text is emitted as <text> elements, which keeps the files small and
editable but relies on the viewer having the fonts. Pass --text-as-path to
convert glyphs to outlines instead, which is larger but renders identically
everywhere - the safer choice for camera-ready submission.

    uv run python scripts/figures_to_svg.py
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pymupdf

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIG_DIR = REPO_ROOT / "manuscript" / "figures"


def convert_pdf(pdf_path: Path, out_dir: Path, text_as_path: bool = False) -> list[Path]:
    """Convert every page of one PDF to SVG; return the paths written."""
    written: list[Path] = []
    with pymupdf.open(pdf_path) as doc:
        multipage = doc.page_count > 1
        for index, page in enumerate(doc, start=1):
            name = (f"{pdf_path.stem}_p{index:02d}.svg" if multipage
                    else f"{pdf_path.stem}.svg")
            out_path = out_dir / name
            out_path.write_text(page.get_svg_image(text_as_path=text_as_path))
            written.append(out_path)
    return written


def convert_all(fig_dir: Path, out_dir: Path, text_as_path: bool = False) -> list[Path]:
    """Convert every PDF in `fig_dir`; return all SVG paths written."""
    pdfs = sorted(fig_dir.glob("*.pdf"))
    if not pdfs:
        raise FileNotFoundError(f"no PDF figures found in {fig_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for pdf_path in pdfs:
        paths = convert_pdf(pdf_path, out_dir, text_as_path)
        total_kb = sum(p.stat().st_size for p in paths) / 1024
        logger.info("%-45s -> %2d SVG (%.0f kB)", pdf_path.name, len(paths), total_kb)
        written.extend(paths)
    return written


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fig-dir", type=Path, default=DEFAULT_FIG_DIR)
    p.add_argument("--out-dir", type=Path, default=None,
                   help="default: <fig-dir>/svg")
    p.add_argument("--text-as-path", action="store_true",
                   help="outline glyphs instead of emitting <text> elements")
    args = p.parse_args()

    out_dir = args.out_dir or args.fig_dir / "svg"
    written = convert_all(args.fig_dir, out_dir, args.text_as_path)
    total_mb = sum(p.stat().st_size for p in written) / 1e6
    logger.info("wrote %d SVG files (%.1f MB) -> %s", len(written), total_mb, out_dir)


if __name__ == "__main__":
    main()
