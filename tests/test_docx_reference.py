"""Cover the styling that decides what the submitted .docx files look like.

pandoc styles a .docx from a reference document and says nothing about the
one it was given, so every failure in this module's territory produces a file
that builds, opens and is wrong: the manuscript in Word's own Aptos with blue
headings looked exactly as finished as the Times one. The patches are anchored
on pandoc's default reference document, which is an upstream file, so the
checks below are mostly about the anchors still holding.
"""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import docx_reference as dr  # noqa: E402

# A reference document's styles.xml in miniature: the two docDefaults blocks
# and one style of each shape the patcher has to handle -- a bare one, one
# with only w:pPr, and one with both property blocks to replace.
STYLES_HEAD = """<?xml version="1.0" encoding="utf-8"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:docDefaults>
    <w:rPrDefault>
      <w:rPr>
        <w:rFonts w:asciiTheme="minorHAnsi" w:hAnsiTheme="minorHAnsi"/>
        <w:sz w:val="24"/>
      </w:rPr>
    </w:rPrDefault>
    <w:pPrDefault>
      <w:pPr>
        <w:spacing w:after="200"/>
      </w:pPr>
    </w:pPrDefault>
  </w:docDefaults>
  <w:style w:type="paragraph" w:styleId="BodyText">
    <w:name w:val="Body Text"/><w:basedOn w:val="Normal"/>
    <w:link w:val="BodyTextChar"/><w:qFormat/>
    <w:pPr><w:spacing w:before="180" w:after="180"/></w:pPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Compact">
    <w:name w:val="Compact"/><w:basedOn w:val="BodyText"/>
    <w:pPr><w:spacing w:before="36" w:after="36"/></w:pPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Title">
    <w:name w:val="Title"/><w:basedOn w:val="Normal"/>
    <w:pPr><w:jc w:val="center"/></w:pPr>
    <w:rPr><w:rFonts w:asciiTheme="majorHAnsi"/><w:sz w:val="56"/></w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Author">
    <w:name w:val="Author"/><w:basedOn w:val="Title"/>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Abstract">
    <w:name w:val="Abstract"/><w:basedOn w:val="Normal"/>
    <w:rPr><w:sz w:val="20"/></w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="AbstractTitle">
    <w:name w:val="Abstract Title"/><w:basedOn w:val="Normal"/>
  </w:style>
"""

HEADINGS = "\n".join(
    f'  <w:style w:type="paragraph" w:styleId="Heading{n}">'
    f'<w:name w:val="heading {n}"/><w:basedOn w:val="Normal"/>'
    f'<w:link w:val="Heading{n}Char"/>'
    f'<w:pPr><w:outlineLvl w:val="{n - 1}"/></w:pPr>'
    f'<w:rPr><w:color w:val="0F4761" w:themeColor="accent1"/></w:rPr>'
    f"</w:style>\n"
    f'  <w:style w:type="character" w:styleId="Heading{n}Char">'
    f'<w:name w:val="Heading {n} Char"/><w:link w:val="Heading{n}"/>'
    f'<w:rPr><w:color w:val="0F4761"/></w:rPr></w:style>'
    for n in range(1, 7)
)

STYLES = STYLES_HEAD + HEADINGS + "\n</w:styles>\n"

THEME = (
    '<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '<a:majorFont><a:latin typeface="Aptos Display"/></a:majorFont>'
    '<a:minorFont><a:latin typeface="Aptos"/></a:minorFont>'
    "</a:theme>"
)


def style_body(xml: str, style_id: str) -> str:
    match = re.search(
        rf'<w:style\b[^>]*w:styleId="{style_id}"[^>]*>(.*?)</w:style>',
        xml, re.DOTALL)
    assert match is not None, f"no style {style_id}"
    return match.group(1)


@pytest.fixture()
def patched() -> str:
    return dr.patch_styles(STYLES)


def test_document_defaults_name_times_rather_than_the_theme(patched: str):
    defaults = re.search(r"<w:rPrDefault>.*?</w:rPrDefault>", patched, re.DOTALL)
    assert defaults is not None
    assert 'w:ascii="Times New Roman"' in defaults.group(0)
    # The theme reference has to go, not just be joined: Word resolves
    # asciiTheme in preference to nothing, and leaving both is how a file
    # comes out half in Aptos.
    assert "minorHAnsi" not in defaults.group(0)


def test_body_text_is_justified_and_not_indented(patched: str):
    body = style_body(patched, "BodyText")
    assert '<w:jc w:val="both"/>' in body
    assert '<w:ind w:firstLine="0"/>' in body


def test_body_text_keeps_a_gap_between_paragraphs(patched: str):
    # With neither an indent nor a gap, consecutive paragraphs run together
    # with nothing to separate them, so the two settings arrive as a pair.
    assert f'w:after="{dr.PARA_SPACE}"' in style_body(patched, "BodyText")
    assert dr.PARA_SPACE > 0


def test_table_cells_are_neither_justified_nor_indented(patched: str):
    # Compact is what pandoc puts in table cells. It is basedOn BodyText, so
    # both properties are inherited unless they are overridden here.
    body = style_body(patched, "Compact")
    assert '<w:jc w:val="left"/>' in body
    assert '<w:ind w:firstLine="0"/>' in body


@pytest.mark.parametrize("level", range(1, 7))
def test_headings_are_twelve_point_bold_and_black(patched: str, level: int):
    body = style_body(patched, f"Heading{level}")
    assert 'w:ascii="Times New Roman"' in body
    assert "<w:b/>" in body
    assert f'<w:sz w:val="{dr.BODY_SZ}"/>' in body
    assert '<w:color w:val="000000"/>' in body
    assert "themeColor" not in body


@pytest.mark.parametrize("level", range(1, 7))
def test_headings_keep_their_outline_level(patched: str, level: int):
    # w:outlineLvl is what Word's navigation pane and any generated table of
    # contents read. Replacing a heading's w:pPr without re-emitting it leaves
    # a document whose headings look right and outline as body text.
    assert f'<w:outlineLvl w:val="{level - 1}"/>' in style_body(
        patched, f"Heading{level}")


def test_the_linked_character_styles_match_their_headings(patched: str):
    for level in range(1, 7):
        char = style_body(patched, f"Heading{level}Char")
        assert 'w:ascii="Times New Roman"' in char
        assert "<w:b/>" in char
        assert '<w:color w:val="0F4761"' not in char


def test_only_the_subsection_levels_are_italic(patched: str):
    # pandoc does not carry LaTeX's section numbering into a .docx, so italic
    # is what keeps a uniformly bold Heading 1 and Heading 2 apart.
    assert "<w:i/>" not in style_body(patched, "Heading1")
    for level in range(2, 7):
        assert "<w:i/>" in style_body(patched, f"Heading{level}")


def test_the_author_line_is_not_bold(patched: str):
    # Author is basedOn Title, so the bold has to be turned off explicitly
    # rather than simply left out.
    assert '<w:b w:val="0"/>' in style_body(patched, "Author")


def test_the_abstract_is_set_at_the_body_size(patched: str):
    body = style_body(patched, "Abstract")
    assert f'<w:sz w:val="{dr.BODY_SZ}"/>' in body
    assert '<w:sz w:val="20"/>' not in body


def test_set_style_keeps_the_metadata_it_does_not_own(patched: str):
    # w:name, w:basedOn and w:link are pandoc's contract with the document it
    # writes; a style rewritten whole would break the inheritance the rest of
    # the sheet is built on.
    body = style_body(patched, "BodyText")
    assert '<w:name w:val="Body Text"/>' in body
    assert '<w:basedOn w:val="Normal"/>' in body
    assert '<w:link w:val="BodyTextChar"/>' in body


def test_set_style_leaves_exactly_one_of_each_property_block(patched: str):
    body = style_body(patched, "Title")
    assert body.count("<w:pPr>") == 1
    assert body.count("<w:rPr>") == 1
    assert '<w:sz w:val="56"/>' not in body


def test_set_style_refuses_a_style_it_cannot_find():
    with pytest.raises(ValueError, match="Absent"):
        dr.set_style(STYLES, "Absent", "<w:pPr/>", "")


def test_patch_styles_refuses_a_reference_without_the_default_fonts():
    # pandoc's default reference document moves between releases. A patch that
    # silently stopped matching would drop the one rule that makes every
    # unnamed style Times, out of a build that still reports success.
    stripped = re.sub(r"<w:rFonts\b[^>]*/>", "", STYLES, count=1)
    with pytest.raises(ValueError, match="default fonts"):
        dr.patch_styles(stripped)


def test_patch_theme_refuses_a_theme_missing_a_face():
    with pytest.raises(ValueError, match="found 1"):
        dr.patch_theme(THEME.replace("<a:minorFont>", "<a:otherFont>"))


def test_patch_theme_points_both_faces_at_times():
    assert re.findall(r'typeface="([^"]*)"', dr.patch_theme(THEME)) == [
        dr.FONT, dr.FONT]


def _pandoc_missing() -> bool:
    try:
        dr._pandoc()
    except (ImportError, OSError):
        return True
    return False


@pytest.mark.skipif(_pandoc_missing(),
                    reason="pandoc is build-time only, not a project dependency")
def test_the_built_reference_is_a_docx_carrying_the_patches(tmp_path: Path):
    out = dr.build_reference(tmp_path / "reference.docx")
    with zipfile.ZipFile(out) as archive:
        names = set(archive.namelist())
        styles = archive.read("word/styles.xml").decode()
        theme = archive.read("word/theme/theme1.xml").decode()
    # The parts pandoc needs have to survive being rewritten, not just the
    # two that were patched.
    assert {"word/document.xml", "word/settings.xml",
            "[Content_Types].xml"} <= names
    assert 'w:ascii="Times New Roman"' in styles
    assert "Aptos" not in theme
