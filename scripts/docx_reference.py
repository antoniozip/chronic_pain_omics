"""Build the pandoc reference document the submitted .docx files are styled from.

pandoc takes every paragraph and character style from a reference document,
and its stock one is Word's own look: the Aptos theme font, a 28 pt title,
headings in accent blue, ragged right. Nothing in the LaTeX preamble reaches
it -- that governs the PDF alone -- so the Word files the journal receives had
never carried the paper's typography, and changing the preamble to fix it
changes nothing here.

The reference is derived from pandoc's own default rather than tracked as a
binary, so every styling decision is readable in the diff instead of being
buried in a zip. That makes the patch depend on pandoc's default, which moves
between releases, so **every substitution asserts that it matched**: a patch
that silently missed would leave the styling absent from a file that still
builds, reports success, and looks finished.

Build-time only, like build_manuscript_docx.py: pandoc is reached through
pypandoc when it is not on PATH.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)

#: Times New Roman at 12 pt throughout. w:sz is in half-points, so 12 pt is 24
#: and 14 pt is 28; lengths are twips (1/20 pt), so the 6 pt between paragraphs
#: is 120.
FONT = "Times New Roman"
BODY_SZ = 24
TITLE_SZ = 28
PARA_SPACE = 120

_RFONTS = (f'<w:rFonts w:ascii="{FONT}" w:hAnsi="{FONT}" w:cs="{FONT}"'
           f' w:eastAsia="{FONT}"/>')
_BLACK = '<w:color w:val="000000"/>'
_SZ = f'<w:sz w:val="{BODY_SZ}"/><w:szCs w:val="{BODY_SZ}"/>'


def _pandoc() -> str:
    """The pandoc executable, from PATH or from the pypandoc wheel."""
    exe = shutil.which("pandoc")
    if exe:
        return exe
    # Build-time only; see the module docstring.
    import pypandoc

    return pypandoc.get_pandoc_path()


def _sub_once(pattern: re.Pattern[str], repl: str, xml: str, what: str) -> str:
    """Apply a substitution that must match exactly once.

    Raising on 0 or 2 matches is the whole point: pandoc's default reference
    document is an upstream file, and a pattern that stopped matching would
    otherwise drop one styling rule out of a build that still succeeds.
    """
    n = len(pattern.findall(xml))
    if n != 1:
        raise ValueError(
            f"{what}: expected exactly one match in pandoc's default reference "
            f"document, found {n}. pandoc's default has probably changed.")
    return pattern.sub(repl, xml, count=1)


def _ppr(*parts: str) -> str:
    """A w:pPr from its children, which must already be in schema order."""
    return "<w:pPr>" + "".join(parts) + "</w:pPr>"


def _rpr(*parts: str) -> str:
    """A w:rPr from its children, which must already be in schema order."""
    return "<w:rPr>" + "".join(parts) + "</w:rPr>"


def _heading(level: int, before: int, after: int, italic: bool) -> tuple[str, str]:
    """One heading level: Times, 12 pt, bold, black, left.

    Every level is the same size and weight, which is what was asked for.
    Levels below the first add italic because pandoc does not carry LaTeX's
    section numbering into a .docx: without it, a uniformly bold Heading 1 and
    Heading 2 are indistinguishable on the page, and the paper runs two levels
    deep throughout the Results and Methods.

    w:outlineLvl is 0-based and must be kept: it is what Word's navigation
    pane and any generated table of contents read.
    """
    ppr = _ppr('<w:keepNext/><w:keepLines/>',
               f'<w:spacing w:before="{before}" w:after="{after}"/>',
               '<w:jc w:val="left"/>',
               f'<w:outlineLvl w:val="{level - 1}"/>')
    rpr = _rpr(_RFONTS, '<w:b/>', '<w:i/>' if italic else '', _BLACK, _SZ)
    return ppr, rpr


def _style_table() -> dict[str, tuple[str, str]]:
    """The w:pPr and w:rPr to force on each style, keyed by w:styleId.

    Only the styles the manuscript and supplementary actually use are listed;
    everything else inherits Times 12 pt from the document defaults.
    """
    styles: dict[str, tuple[str, str]] = {
        # Justified, no first-line indent, and 6 pt between paragraphs.
        # The space is not decoration: with neither an indent nor a gap,
        # consecutive paragraphs run together with nothing to separate them,
        # so the two settings come as a pair.
        "BodyText": (
            _ppr(f'<w:spacing w:before="0" w:after="{PARA_SPACE}" w:line="240"'
                 ' w:lineRule="auto"/>',
                 '<w:ind w:firstLine="0"/>',
                 '<w:jc w:val="both"/>'),
            ""),
        # Compact is what pandoc puts in table cells -- 224 of them here --
        # where an inherited first-line indent and justification both read as
        # damage.
        "Compact": (
            _ppr('<w:spacing w:before="36" w:after="36"/>',
                 '<w:ind w:firstLine="0"/>',
                 '<w:jc w:val="left"/>'),
            ""),
        "Title": (
            _ppr('<w:keepNext/><w:keepLines/>',
                 '<w:spacing w:before="0" w:after="240" w:line="240"'
                 ' w:lineRule="auto"/>',
                 '<w:contextualSpacing/>',
                 '<w:jc w:val="center"/>'),
            _rpr(_RFONTS, '<w:b/>', _BLACK,
                 f'<w:sz w:val="{TITLE_SZ}"/><w:szCs w:val="{TITLE_SZ}"/>')),
        # Author is basedOn Title, so the bold has to be turned off rather
        # than left out.
        "Author": (
            _ppr('<w:keepNext/><w:keepLines/>',
                 '<w:spacing w:before="0" w:after="120"/>',
                 '<w:jc w:val="center"/>'),
            _rpr('<w:b w:val="0"/>', _SZ)),
        "AbstractTitle": (
            _ppr('<w:keepNext/><w:keepLines/>',
                 '<w:spacing w:before="240" w:after="60"/>',
                 '<w:jc w:val="left"/>'),
            _rpr(_RFONTS, '<w:b/>', _BLACK, _SZ)),
        # pandoc sets the abstract two points smaller than the body; the
        # request was one size throughout.
        "Abstract": (
            _ppr('<w:spacing w:before="0" w:after="120"/>',
                 '<w:ind w:firstLine="0"/>',
                 '<w:jc w:val="both"/>'),
            _rpr(_SZ)),
    }
    spacing = {1: (240, 120), 2: (200, 80), 3: (200, 80),
               4: (160, 40), 5: (160, 40), 6: (160, 40)}
    for level, (before, after) in spacing.items():
        ppr, rpr = _heading(level, before, after, italic=level > 1)
        styles[f"Heading{level}"] = (ppr, rpr)
        # The linked character style carries the same run properties, and is
        # what Word applies when a heading is selected as text.
        styles[f"Heading{level}Char"] = ("", rpr)
    return styles


# A style's own property blocks, self-closing or not. \b keeps this off
# w:pPrDefault and w:rPrDefault, which are children of w:docDefaults and are
# never inside a w:style.
_PR_BLOCK_RE = re.compile(
    r"<w:pPr\b(?:[^>]*/>|.*?</w:pPr>)|<w:rPr\b(?:[^>]*/>|.*?</w:rPr>)",
    re.DOTALL)


def set_style(xml: str, style_id: str, ppr: str, rpr: str) -> str:
    """Replace one style's w:pPr and w:rPr, keeping the rest of it.

    w:name, w:basedOn, w:next and w:link are pandoc's contract with the
    document it writes -- rewriting a style whole would break the inheritance
    the other styles rely on -- so only the two property blocks are touched.
    Schema order puts both after all the metadata, which is why the new blocks
    are appended rather than substituted where the old ones were.
    """
    pattern = re.compile(
        r'(<w:style\b[^>]*w:styleId="' + re.escape(style_id) + r'"[^>]*>)'
        r"(.*?)(</w:style>)", re.DOTALL)

    def repl(match: re.Match[str]) -> str:
        body = _PR_BLOCK_RE.sub("", match.group(2)).rstrip()
        return match.group(1) + body + ppr + rpr + match.group(3)

    return _sub_once(pattern, repl, xml, f"style {style_id}")


_DEFAULT_FONTS_RE = re.compile(
    r"(<w:rPrDefault>\s*<w:rPr>\s*)<w:rFonts\b[^>]*/>")
_DEFAULT_SPACING_RE = re.compile(
    r"(<w:pPrDefault>\s*<w:pPr>\s*)<w:spacing\b[^>]*/>")


def patch_styles(xml: str) -> str:
    """Times New Roman throughout, justified body, 12 pt bold headings."""
    # The document defaults are what make every style Times, including the
    # ones below that are never named again.
    xml = _sub_once(_DEFAULT_FONTS_RE, r"\1" + _RFONTS, xml, "default fonts")
    xml = _sub_once(_DEFAULT_SPACING_RE,
                    r'\1<w:spacing w:after="0" w:line="240" w:lineRule="auto"/>',
                    xml, "default paragraph spacing")
    for style_id, (ppr, rpr) in _style_table().items():
        xml = set_style(xml, style_id, ppr, rpr)
    return xml


_THEME_LATIN_RE = re.compile(r'(<a:(?:major|minor)Font>\s*<a:latin) typeface="[^"]*"')


def patch_theme(xml: str) -> str:
    """Point the theme's major and minor faces at Times as well.

    Redundant against the explicit w:rFonts above, and kept because a style
    this module does not name can still reach the theme, and Word resolves
    asciiTheme before it falls back to anything.
    """
    n = len(_THEME_LATIN_RE.findall(xml))
    if n != 2:
        raise ValueError(f"expected a major and a minor latin face, found {n}")
    return _THEME_LATIN_RE.sub(rf'\1 typeface="{FONT}"', xml)


PATCHES = {"word/styles.xml": patch_styles, "word/theme/theme1.xml": patch_theme}


def build_reference(dest: Path) -> Path:
    """Write the patched reference document, returning its path."""
    default = subprocess.run(
        [_pandoc(), "--print-default-data-file", "reference.docx"],
        capture_output=True, check=True).stdout
    source = dest.with_name(dest.stem + "_pandoc_default.docx")
    source.write_bytes(default)

    applied = set()
    with zipfile.ZipFile(source) as zin, \
            zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            patch = PATCHES.get(item.filename)
            if patch is not None:
                data = patch(data.decode("utf-8")).encode("utf-8")
                applied.add(item.filename)
            zout.writestr(item, data)
    missing = set(PATCHES) - applied
    if missing:
        raise ValueError(f"pandoc's default reference document has no {missing}")
    logger.info("built reference document %s (%d bytes)", dest.name,
                dest.stat().st_size)
    return dest
