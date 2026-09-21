"""Render manuscript/manuscript.tex to Word (.docx).

pandoc cannot embed PDF figures in a docx, so each figure is first rasterized
to PNG with PyMuPDF and the LaTeX source is rewritten to point at the PNGs in a
scratch copy. The manuscript source itself is never modified.

Build-time only: pypandoc-binary and pymupdf are not project dependencies.

    uv run python scripts/build_manuscript_docx.py
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import tempfile
from pathlib import Path

import docx_reference

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TEX = REPO_ROOT / "manuscript" / "manuscript.tex"
RASTER_DPI = 200

_INCLUDEGRAPHICS_RE = re.compile(r"(\\includegraphics(?:\[[^\]]*\])?\{)([^}]+)(\})")
# The abstract sits inside \twocolumn[...], which pandoc discards wholesale,
# so the docx was published with no abstract at all. Unwrap the block and
# rename onecolabstract to the standard environment pandoc understands.
_TWOCOLUMN_RE = re.compile(r"\\twocolumn\[(.*?)\n\]", re.DOTALL)
# `\captionof` comes from the caption package and is used for the graphical
# abstract, which sits inside \twocolumn[...] where a float is not allowed.
# Pandoc does not know the macro and drops its argument silently, so the
# figure arrived in the docx with no caption at all -- the same class of loss
# as the missing abstract above. Rewriting it to a plain braced group keeps
# the text, since the docx has no float numbering to attach it to anyway.
#
# The blank line is load-bearing. Pandoc's LaTeX reader drops the caption
# entirely if it follows \includegraphics without a paragraph break, which is
# how a first attempt at this fix still produced a captionless figure.
_CAPTIONOF_RE = re.compile(r"\\captionof\{figure\}\{", re.DOTALL)
# The references moved from a hand-written thebibliography to BibTeX, and
# pandoc resolves neither \bibliographystyle nor \bibliography: run against
# the raw .tex it produces a docx with no reference list at all, silently.
# BibTeX has already formatted the list into manuscript.bbl, which *is* a
# thebibliography environment, so splicing that in gives the docx the same
# numbered list the PDF carries rather than a second rendering that can differ.
_BIBLIOGRAPHY_RE = re.compile(
    r"\\bibliographystyle\{[^}]*\}\s*\n\\bibliography\{[^}]*\}")


def inline_bbl(source: str, bbl_path: Path, numbers: dict[str, int]) -> str:
    """Replace the BibTeX commands with the list BibTeX produced.

    Raises when the .bbl is missing rather than dropping the references: a
    docx that silently ships without a bibliography is the failure this
    function exists to prevent.
    """
    if not _BIBLIOGRAPHY_RE.search(source):
        return source                      # a source that embeds its own list
    if not bbl_path.exists():
        raise FileNotFoundError(
            f"{bbl_path} not found. Build the bibliography first:\n"
            f"  cd {bbl_path.parent} && pdflatex manuscript && bibtex manuscript"
        )
    # thebibliography prints its own "References" heading in the article class,
    # but pandoc renders the environment as a bare list, so the docx carried an
    # unheaded block of 59 numbered entries. Supply the heading the PDF shows.
    body = ("\\section*{References}\n\n"
            + bbl_as_numbered_paragraphs(bbl_path.read_text(), numbers))
    return _BIBLIOGRAPHY_RE.sub(lambda _: body, source)


# Pandoc does not implement `thebibliography`: it drops the \begin and leaves
# the widest-label argument behind as a stray paragraph ("10"), and it spills
# elsarticle-num's \expandafter/\csname guards into the text as the words
# "url urlprefix href".
#
# An enumerate fixes that but lets Word number the list itself, and a Word
# list renumbers from whatever the document's list state is rather than from
# the \bibcite map the in-text markers come from. Each entry is therefore a
# plain paragraph carrying its printed label as text: "1.", the AMA form the
# PDF prints since the move to Brain Communications (it was "[1]" under
# elsarticle-num).
#
# The guards cannot simply be deleted: they define \path, and \path is what
# carries the DOI inside each \href. Removing them took all 59 DOIs out of the
# docx while leaving the reference list looking complete. They are restated as
# \providecommand, which pandoc does expand.
_BBL_HEAD_RE = re.compile(
    r"\A\\begin\{thebibliography\}\{[^}]*\}.*?(?=\\bibitem)", re.DOTALL)
_BBL_PROVIDES = (
    "\\providecommand{\\url}[1]{\\texttt{#1}}\n"
    "\\providecommand{\\urlprefix}{URL }\n"
    "\\providecommand{\\path}[1]{#1}\n"
)
_BIBITEM_RE = re.compile(r"\\bibitem\{([^}]*)\}")


def bbl_as_numbered_paragraphs(bbl: str, numbers: dict[str, int]) -> str:
    """Recast a BibTeX .bbl as paragraphs labelled the way the PDF labels them.

    Labels come from the same \bibcite map the in-text markers are built from,
    so the two cannot disagree; a key the .aux does not number would mean the
    build is a pass behind, and is refused rather than silently unlabelled.
    """
    # Lambda replacements: the replacement text is full of backslashes, and
    # re.sub reads those as group escapes ("bad escape \p" for \providecommand).
    body = _BBL_HEAD_RE.sub(lambda _: _BBL_PROVIDES, bbl)

    def label(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in numbers:
            raise KeyError(
                f"{key} is in the .bbl but carries no \\bibcite number; "
                "re-run pdflatex after bibtex"
            )
        return f"{numbers[key]}.~"

    body = _BIBITEM_RE.sub(label, body)
    # Pandoc treats \newblock as a command with an argument and drops the group
    # that follows, so a title opening with a protected word ("{NCBI GEO}:
    # archive ...") lost it. The blocks are only spacing; a space does as well.
    body = re.sub(r"\\newblock\s*", " ", body)
    return body.replace("\\end{thebibliography}", "")


# Pandoc parses \cite into a citation element and then, with no bibliography
# metadata and no --citeproc, renders it as nothing at all: the docx carried
# the reference list but not one marker pointing into it, and had done so
# since before the move to BibTeX. Substituting the numbers ourselves keeps
# the docx and the PDF on one source of truth, since \bibcite in the .aux is
# what LaTeX numbered the printed list from.
_CITE_RE = re.compile(r"\\cite\{([^}]*)\}")
# A citation with the space before it and the one punctuation mark after it
# that the journal places *before* the number. Brain Communications sets the
# superscript after a full stop or comma, before a semicolon or colon, with no
# space in front -- which is what cite.sty prints with \CiteMoveChars set to
# `.,`. The docx has to do the same, or the two artifacts disagree about where
# every citation sits.
# The whitespace may include the one newline a wrapped source line leaves
# before a \cite; a blank line (a paragraph break) is never consumed.
_CITE_IN_CONTEXT_RE = re.compile(r"[ \t]*\n?[ \t]*\\cite\{([^}]*)\}([.,]?)")
_BIBCITE_RE = re.compile(r"\\bibcite\{([^}]*)\}\{\{?(\d+)")


def citation_numbers(aux_path: Path) -> dict[str, int]:
    """Key -> printed number, read from LaTeX's own \\bibcite records."""
    if not aux_path.exists():
        raise FileNotFoundError(
            f"{aux_path} not found. Build the bibliography first:\n"
            f"  cd {aux_path.parent} && pdflatex manuscript && bibtex manuscript"
            " && pdflatex manuscript"
        )
    return {k: int(n) for k, n in _BIBCITE_RE.findall(aux_path.read_text())}


def _compress(numbers: list[int]) -> str:
    """Render a sorted citation list the way the journal writes it: 1,3,6-9."""
    out: list[str] = []
    i = 0
    while i < len(numbers):
        j = i
        while j + 1 < len(numbers) and numbers[j + 1] == numbers[j] + 1:
            j += 1
        # A run of two prints as two entries; cite.sty only bridges three or more.
        if j - i >= 2:
            out.append(f"{numbers[i]}-{numbers[j]}")
        else:
            out.extend(str(n) for n in numbers[i:j + 1])
        i = j + 1
    # No space after the comma: the journal's own example is "1,3,6-9", and
    # \citepunct is set to match in the PDF.
    return ",".join(out)


def numeric_citations(source: str, numbers: dict[str, int]) -> str:
    """Replace every \\cite with the superscript number the PDF prints."""
    if not _CITE_RE.search(source):
        return source

    def render(match: re.Match[str]) -> str:
        keys = [k.strip() for k in match.group(1).split(",") if k.strip()]
        missing = [k for k in keys if k not in numbers]
        if missing:
            raise KeyError(
                f"{', '.join(missing)} is cited but carries no \\bibcite "
                "number; re-run pdflatex after bibtex"
            )
        marker = _compress(sorted({numbers[k] for k in keys}))
        return match.group(2) + "\\textsuperscript{" + marker + "}"

    text, n = _CITE_IN_CONTEXT_RE.subn(render, source)
    logger.info("rendered %d citation(s) as numbers", n)
    return text


# xr resolves \\ref against another document's .aux at LaTeX time, which pandoc
# cannot do: the supplementary .docx rendered "Figure [fig:forest]" where the
# PDF prints "Figure 3". Same remedy as the citations -- read the numbers
# LaTeX recorded and substitute them.
# The optional [prefix] is xr's own: every label imported from that document
# is known here as prefix + label, which is how the supplementary keeps the
# paper's \bibcite entries from colliding with its own reference list.
_EXTERNALDOC_RE = re.compile(r"\\externaldocument(?:\[([^\]]*)\])?\{([^}]+)\}")
_NEWLABEL_RE = re.compile(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}")
_LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
_REF_RE = re.compile(r"\\ref\{([^}]+)\}")


def label_numbers(aux_path: Path) -> dict[str, str]:
    """Label -> printed number, from LaTeX's own \\newlabel records."""
    if not aux_path.exists():
        raise FileNotFoundError(
            f"{aux_path} not found; build the document it belongs to first"
        )
    return dict(_NEWLABEL_RE.findall(aux_path.read_text()))


def resolve_external_refs(source: str, base_dir: Path) -> str:
    """Substitute \\ref numbers for labels defined in another document.

    Only labels the source does not define itself are touched; pandoc resolves
    its own \\label/\\ref pairs correctly, and rewriting those would replace a
    live cross-reference with frozen text.
    """
    external = _EXTERNALDOC_RE.findall(source)
    if not external:
        return source
    own = set(_LABEL_RE.findall(source))
    numbers: dict[str, str] = {}
    for prefix, name in external:
        numbers.update({prefix + key: value for key, value in
                        label_numbers((base_dir / name).with_suffix(".aux")).items()})

    def render(match: re.Match[str]) -> str:
        key = match.group(1)
        if key in own:
            return match.group(0)
        if key not in numbers:
            raise KeyError(
                f"{key} is referenced but appears in no \\externaldocument's "
                ".aux; rebuild that document first"
            )
        return numbers[key]

    text, n = _REF_RE.subn(render, source)
    logger.info("resolved %d external cross-reference(s)", n)
    return text


# Pandoc attaches a figure caption only when \caption *follows* the graphic.
# Every figure caption in this manuscript sits above its graphic, which is what
# the journal asks for, and read literally that drops all five captions from
# the .docx without a word -- the same silent loss as the \captionof case
# above. The caption is therefore lifted out of the float and left as an
# ordinary paragraph in front of it, which is the position the source asks for
# and the one arrangement pandoc cannot reinterpret.
#
# Both anchors are line-start: the preamble carries the words
# "\begin{figure*} for the volcano comparison" inside a comment, and without
# ^ the scan starts there and swallows four figures into one match.
_FIGURE_ENV_RE = re.compile(
    r"^\\begin\{(figure\*?)\}(.*?)^\\end\{\1\}\n", re.DOTALL | re.MULTILINE)
_CAPTION_OPEN = "\\caption{"


def _brace_end(text: str, open_brace: int) -> int:
    """Index just past the brace group that opens at `open_brace`."""
    depth, i = 1, open_brace + 1
    while depth:
        if text[i] == "\\":                  # an escaped brace is not a brace
            i += 2
            continue
        depth += (text[i] == "{") - (text[i] == "}")
        i += 1
    return i


def hoist_figure_captions(source: str) -> str:
    """Lift a caption that precedes its graphic out in front of the float.

    Only captions that come first are touched. One that follows the graphic is
    left alone, since that is the arrangement pandoc already reads correctly.
    """
    hoisted = 0

    def rewrite(match: re.Match[str]) -> str:
        nonlocal hoisted
        kind, body = match.group(1), match.group(2)
        cap = body.find(_CAPTION_OPEN)
        img = body.find(r"\includegraphics")
        if cap < 0 or (0 <= img < cap):
            return match.group(0)
        end = _brace_end(body, cap + len(_CAPTION_OPEN) - 1)
        text = body[cap + len(_CAPTION_OPEN):end - 1]
        hoisted += 1
        return (f"{text}\n\n"
                f"\\begin{{{kind}}}{body[:cap] + body[end:]}\\end{{{kind}}}\n")

    text = _FIGURE_ENV_RE.sub(rewrite, source)
    logger.info("hoisted %d figure caption(s) pandoc would have dropped", hoisted)
    return text


# \path, and the \repopath wrapper around it, never reach a Word file: pandoc
# drops them, and the supplementary shipped "reason codes ()" where the PDF
# names literature/prisma/. As \texttt they survive, with underscores escaped
# because \path took its argument verbatim.
_REPOPATH_RE = re.compile(r"\\(?:repopath|path)\{([^{}]*)\}")


def expand_repository_paths(source: str) -> str:
    """Rewrite each \\repopath{...} or \\path{...} as \\texttt{...}."""
    def code(match: re.Match[str]) -> str:
        return "\\texttt{" + match.group(1).replace("_", "\\_") + "}"

    return _REPOPATH_RE.sub(code, source)


_NUMBERED_CAPTION_RE = re.compile(r"\\caption\{")


def numbers_its_captions(source: str) -> bool:
    """Whether the source has captions for pandoc to number natively.

    The supplementary labels its one printed table itself (\\caption*), and
    native numbering would call it "Table 1".
    """
    return bool(_NUMBERED_CAPTION_RE.search(source))


def unwrap_abstract(source: str) -> str:
    """Make the abstract an ordinary headed section, where the source puts it.

    Two things would otherwise lose or misplace it. Inside \\twocolumn[...]
    (the Journal of Pain layout) pandoc discards it wholesale. As an abstract
    environment pandoc lifts it into metadata and prints it straight after the
    author line, ahead of the affiliations, running title and correspondence
    that Brain Communications wants on the title page.
    """
    def unwrap(match: re.Match[str]) -> str:
        return match.group(1)

    text, _ = _TWOCOLUMN_RE.subn(unwrap, source)
    text = text.replace(r"\begin{onecolabstract}", r"\begin{abstract}")
    text = text.replace(r"\end{onecolabstract}", r"\end{abstract}")
    if r"\begin{abstract}" not in text:
        # Expected for the supplementary document, which has none; for the
        # paper it means the abstract has gone missing from the Word file.
        logger.warning("no abstract environment found; the docx has no abstract "
                       "(expected only for the supplementary document)")
    text = re.sub(r"[ \t]*\\begin\{abstract\}[ \t]*\n?", lambda _: "\\section*{Abstract}\n\n", text)
    text = re.sub(r"[ \t]*\\end\{abstract\}[ \t]*", "", text)
    text, n_cap = _CAPTIONOF_RE.subn("\n\n{", text)
    if n_cap:
        logger.info("rewrote %d \\captionof block(s) pandoc would have dropped", n_cap)
    return text


def rasterize(pdf_path: Path, out_path: Path, dpi: int = RASTER_DPI) -> Path:
    """Render page 1 of a PDF figure to PNG."""
    # Build-time only; see the module docstring.
    import pymupdf

    with pymupdf.open(pdf_path) as doc:
        page = doc.load_page(0)
        page.get_pixmap(dpi=dpi).save(out_path)
    return out_path


# Pandoc drops an \\input it cannot resolve, warning as it goes and producing a
# document without that section: the first supplementary .docx came out empty
# for exactly this reason, and copying the file next to the staged source did
# not fix it. Substitute the text in ourselves, as the .bbl already is, so the
# staged source is self-contained and nothing depends on how pandoc searches.
_INPUT_RE = re.compile(r"^[ \t]*\\(?:input|include)\{([^}]+)\}[ \t]*$",
                       re.MULTILINE)


def inline_includes(source: str, base_dir: Path, _depth: int = 0) -> str:
    """Replace each \\input line with the file's text, recursively."""
    if _depth > 8:
        raise RecursionError("\\input nested more than 8 deep; is there a cycle?")

    def substitute(match: re.Match[str]) -> str:
        path = base_dir / match.group(1)
        if path.suffix != ".tex":
            path = path.with_suffix(".tex")
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is \\input by the source but does not exist; the docx "
                "would be missing that section"
            )
        logger.info("inlined %s", path.name)
        return inline_includes(path.read_text(), base_dir, _depth + 1)

    return _INPUT_RE.sub(substitute, source)


def stage_sources(tex_path: Path, work_dir: Path) -> Path:
    """Copy the .tex into work_dir, PNG figures and BibTeX list inlined."""
    source = inline_includes(tex_path.read_text(), tex_path.parent)
    source = resolve_external_refs(source, tex_path.parent)
    fig_dir = work_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    def swap(match: re.Match[str]) -> str:
        prefix, ref, suffix = match.groups()
        original = (tex_path.parent / ref).resolve()
        candidate = original if original.suffix else original.with_suffix(".pdf")
        if not candidate.exists():
            logger.warning("figure not found, leaving reference as-is: %s", ref)
            return match.group(0)
        png = fig_dir / (candidate.stem + ".png")
        rasterize(candidate, png)
        logger.info("rasterized %s -> %s", candidate.name, png.name)
        return f"{prefix}figures/{png.name}{suffix}"

    staged = work_dir / tex_path.name
    # One \bibcite map for both the list labels and the in-text markers, so a
    # reader following "[8]" into the references cannot land on a different
    # entry than the PDF sends them to. Read only when the source needs it:
    # the supplementary document cites nothing and has no bibliography, and
    # demanding an .aux from it would fail the build for no reason.
    needs_numbers = bool(_CITE_RE.search(source)
                         or _BIBLIOGRAPHY_RE.search(source))
    numbers = citation_numbers(tex_path.with_suffix(".aux")) if needs_numbers else {}
    body = inline_bbl(source, tex_path.with_suffix(".bbl"), numbers)
    body = numeric_citations(body, numbers)
    body = expand_repository_paths(body)
    staged.write_text(
        hoist_figure_captions(
            unwrap_abstract(_INCLUDEGRAPHICS_RE.sub(swap, body))))
    return staged


def build(tex_path: Path, out_path: Path) -> Path:
    # Build-time only; see the module docstring.
    import pypandoc

    with tempfile.TemporaryDirectory() as tmp:
        work_dir = Path(tmp)
        staged = stage_sources(tex_path, work_dir)
        # Without --reference-doc the file comes out in Word's own look --
        # Aptos, blue headings, ragged right -- because nothing in the LaTeX
        # preamble reaches a .docx. The reference is built here rather than
        # tracked so that the styling is readable as source.
        reference = docx_reference.build_reference(work_dir / "reference.docx")
        # native_numbering labels each caption "Table 1:" / "Figure 1:" with a
        # Word field. Without it the captions carried no number at all while
        # the text cited "Table 1" and "Fig. 2" -- true of the Journal of Pain
        # submission too, found on the move to Brain Communications.
        numbered = numbers_its_captions(staged.read_text())
        pypandoc.convert_file(
            str(staged), "docx+native_numbering" if numbered else "docx",
            format="latex", outputfile=str(out_path),
            extra_args=["--resource-path", str(work_dir),
                        "--reference-doc", str(reference)],
        )
    logger.info("wrote %s (%.1f MB)", out_path, out_path.stat().st_size / 1e6)
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()
    out = args.out or args.tex.with_suffix(".docx")
    import pypandoc

    if shutil.which("pandoc") is None and not hasattr(pypandoc, "convert_file"):
        raise RuntimeError("pandoc unavailable; install pypandoc-binary")
    build(args.tex, out)


if __name__ == "__main__":
    main()
