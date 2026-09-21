"""Cover the text transforms that decide what reaches the submitted .docx.

Every one of these guards a failure that produced a plausible-looking file.
Pandoc does not implement `thebibliography` or resolve `\\cite`, and it says
nothing when it drops them: the docx shipped with no in-text citation markers
at all until 2026-09-07, and a first attempt at fixing the reference list took
all 59 DOIs out of it while leaving the entries in place. Neither showed up in
the build log, so the checks live here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import build_manuscript_docx as bmd  # noqa: E402

AUX = """\\relax
\\citation{alpha}
\\bibcite{alpha}{{1}{}{{}}{{}}}
\\bibcite{beta}{{2}{}{{}}{{}}}
\\bibcite{gamma}{{3}{}{{}}{{}}}
\\bibcite{delta}{{4}{}{{}}{{}}}
\\bibcite{omega}{{9}{}{{}}{{}}}
"""

BBL = """\\begin{thebibliography}{10}
\\expandafter\\ifx\\csname url\\endcsname\\relax
  \\def\\url#1{\\texttt{#1}}\\fi
\\expandafter\\ifx\\csname urlprefix\\endcsname\\relax\\def\\urlprefix{URL }\\fi

\\bibitem{alpha}
A.~Author, A title, A Journal 1 (2020) 1--9.
\\newblock \\href {https://doi.org/10.1000/a} {\\path{doi:10.1000/a}}.

\\bibitem{beta}
B.~Author, Another title, B Journal 2 (2021) 10--19.
\\newblock \\href {https://doi.org/10.1000/b} {\\path{doi:10.1000/b}}.

\\end{thebibliography}
"""


@pytest.fixture()
def aux(tmp_path: Path) -> Path:
    path = tmp_path / "manuscript.aux"
    path.write_text(AUX)
    return path


def test_citation_numbers_reads_the_aux(aux: Path):
    assert bmd.citation_numbers(aux) == {
        "alpha": 1, "beta": 2, "gamma": 3, "delta": 4, "omega": 9,
    }


def test_citation_numbers_demands_the_aux(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="bibtex"):
        bmd.citation_numbers(tmp_path / "absent.aux")


@pytest.mark.parametrize(("numbers", "expected"), [
    ([1], "1"),
    ([1, 2], "1,2"),                     # cite.sty bridges three, not two
    ([1, 2, 3], "1-3"),
    ([1, 2, 3, 7], "1-3,7"),
    ([1, 3, 5], "1,3,5"),
    ([4, 5, 6, 9, 10], "4-6,9,10"),
])
def test_compress_matches_the_journal_style(numbers, expected):
    """Brain Communications writes a citation list as `1,3,6-9`."""
    assert bmd._compress(numbers) == expected


def _sup(text: str) -> str:
    return "\\textsuperscript{" + text + "}"


def test_numeric_citations_renders_the_printed_marker(aux: Path):
    out = bmd.numeric_citations(r"text \cite{alpha} more \cite{beta, alpha}",
                                bmd.citation_numbers(aux))
    assert out == f"text{_sup('1')} more{_sup('1,2')}"


def test_numeric_citations_sorts_and_bridges_a_run(aux: Path):
    assert bmd.numeric_citations(r"\cite{gamma, alpha, beta}",
                                 bmd.citation_numbers(aux)) == _sup("1-3")


def test_numeric_citations_spans_a_line_break(aux: Path):
    """The manuscript wraps long \\cite lists across lines."""
    assert bmd.numeric_citations("\\cite{alpha,\nomega}",
                                 bmd.citation_numbers(aux)) == _sup("1,9")


@pytest.mark.parametrize(("source", "expected"), [
    ("pain \\cite{alpha}.", "pain." + _sup("1")),
    ("pain \\cite{alpha}, and", "pain," + _sup("1") + " and"),
    ("pain \\cite{alpha}; and", "pain" + _sup("1") + "; and"),
    ("pain \\cite{alpha}: and", "pain" + _sup("1") + ": and"),
    ("(see \\cite{alpha})", "(see" + _sup("1") + ")"),
    # The source wraps lines, so a \\cite often opens one. The newline before
    # it is a space to pandoc, and "population ^1^" is what the first build
    # printed.
    ("population\n\\cite{alpha} and", "population" + _sup("1") + " and"),
    ("PRIDE\n\\cite{alpha}, and", "PRIDE," + _sup("1") + " and"),
])
def test_numeric_citations_follow_the_journal_placement(aux: Path, source, expected):
    """After a full stop or comma, before a semicolon or colon, no space.

    That is the journal's rule and what cite.sty prints in the PDF with
    \\CiteMoveChars set to `.,`; the docx has to agree with it or the two
    artifacts carry the citations in different places.
    """
    assert bmd.numeric_citations(source, bmd.citation_numbers(aux)) == expected


def test_numeric_citations_refuses_an_unnumbered_key(aux: Path):
    """A key absent from the .aux means the build is a pass behind.

    Rendering it as an empty marker would look like prose that never had a
    citation, which is the failure this whole module exists to prevent.
    """
    with pytest.raises(KeyError, match="unknownkey"):
        bmd.numeric_citations(r"\cite{unknownkey}", bmd.citation_numbers(aux))


def test_numeric_citations_leaves_a_source_without_citations_alone():
    """No \\cite means the numbering is never consulted."""
    assert bmd.numeric_citations("plain text", {}) == "plain text"


NUMBERS = {"alpha": 1, "beta": 2}


def test_bbl_keeps_the_doi_macro():
    """\\path carries the DOI inside each \\href.

    Dropping the guard block outright is what removed every DOI from the docx
    while the reference list still looked complete.
    """
    out = bmd.bbl_as_numbered_paragraphs(BBL, NUMBERS)
    assert r"\providecommand{\path}[1]{#1}" in out
    assert out.count(r"\path{doi:") == 2


def test_bbl_drops_what_pandoc_would_spill():
    out = bmd.bbl_as_numbered_paragraphs(BBL, NUMBERS)
    for leak in (r"\expandafter", r"\csname", r"\bibitem", "thebibliography"):
        assert leak not in out, f"{leak} would reach the docx as text"


def test_bbl_labels_entries_the_way_the_pdf_does():
    """The PDF numbers its list "1.", as AMA style does.

    Each entry is a paragraph carrying its label as text, not a Word list, so
    the label cannot be renumbered by Word or drift from the \\bibcite map.
    """
    out = bmd.bbl_as_numbered_paragraphs(BBL, NUMBERS)
    assert "\n1.~" in out and "\n2.~" in out
    assert "{[}" not in out


def test_bbl_refuses_an_entry_the_aux_does_not_number():
    """An unnumbered entry means the build is a pass behind."""
    with pytest.raises(KeyError, match="beta"):
        bmd.bbl_as_numbered_paragraphs(BBL, {"alpha": 1})


def test_inline_bbl_demands_the_bbl(tmp_path: Path):
    """A missing .bbl must stop the build, not silently drop the references."""
    source = "body\n\\bibliographystyle{ama-brain}\n\\bibliography{references}\n"
    with pytest.raises(FileNotFoundError, match="bibtex"):
        bmd.inline_bbl(source, tmp_path / "absent.bbl", NUMBERS)


def test_inline_bbl_supplies_the_heading(tmp_path: Path):
    """thebibliography prints "References" itself; pandoc does not."""
    bbl = tmp_path / "manuscript.bbl"
    bbl.write_text(BBL)
    source = "body\n\\bibliographystyle{ama-brain}\n\\bibliography{references}\n"
    out = bmd.inline_bbl(source, bbl, NUMBERS)
    assert r"\section*{References}" in out
    assert r"\bibliography{references}" not in out


def test_inline_bbl_passes_through_an_embedded_bibliography(tmp_path: Path):
    """A source that already carries its own list needs no .bbl."""
    source = "body\n\\begin{thebibliography}{9}\n\\bibitem{a}A.\n\\end{thebibliography}\n"
    assert bmd.inline_bbl(source, tmp_path / "absent.bbl", NUMBERS) == source


# ---------------------------------------------------------------------------
# \input and external cross-references
# ---------------------------------------------------------------------------
# Both of these fail by producing a document that looks finished. Pandoc drops
# an \input it cannot resolve with only a warning -- the first supplementary
# .docx came out empty that way -- and renders an xr-resolved \ref as the raw
# label, "Figure [fig:forest]" where the PDF prints "Figure 3".


def test_inline_includes_substitutes_the_file(tmp_path: Path):
    (tmp_path / "body.tex").write_text("the body text\n")
    out = bmd.inline_includes("before\n\\input{body}\nafter\n", tmp_path)
    assert "the body text" in out
    assert "\\input" not in out


def test_inline_includes_recurses(tmp_path: Path):
    (tmp_path / "outer.tex").write_text("outer\n\\input{inner}\n")
    (tmp_path / "inner.tex").write_text("inner text\n")
    out = bmd.inline_includes("\\input{outer}\n", tmp_path)
    assert "outer" in out and "inner text" in out


def test_inline_includes_refuses_a_missing_file(tmp_path: Path):
    """Dropping it would ship a docx missing a whole section, silently."""
    with pytest.raises(FileNotFoundError, match="absent"):
        bmd.inline_includes("\\input{absent}\n", tmp_path)


def test_inline_includes_leaves_inputenc_alone(tmp_path: Path):
    """\\usepackage[utf8]{inputenc} contains "input" but is not an include."""
    src = "\\usepackage[utf8]{inputenc}\n"
    assert bmd.inline_includes(src, tmp_path) == src


def _external(tmp_path: Path) -> Path:
    (tmp_path / "main.aux").write_text(
        "\\newlabel{fig:forest}{{3}{5}{}{}{}}\n"
        "\\newlabel{tab:strat}{{1}{4}{}{}{}}\n")
    return tmp_path


def test_resolve_external_refs_substitutes_the_number(tmp_path: Path):
    _external(tmp_path)
    src = "\\externaldocument{main}\nas Figure~\\ref{fig:forest}.\n"
    assert "as Figure~3." in bmd.resolve_external_refs(src, tmp_path)


def test_resolve_external_refs_leaves_own_labels_alone(tmp_path: Path):
    """Pandoc resolves a document's own \\label/\\ref pairs correctly."""
    _external(tmp_path)
    src = ("\\externaldocument{main}\n\\label{fig:mine}\n"
           "see \\ref{fig:mine} and \\ref{fig:forest}\n")
    out = bmd.resolve_external_refs(src, tmp_path)
    assert "\\ref{fig:mine}" in out
    assert "and 3" in out


def test_resolve_external_refs_is_a_no_op_without_externaldocument(tmp_path: Path):
    src = "see \\ref{fig:local}\n"
    assert bmd.resolve_external_refs(src, tmp_path) == src


def test_resolve_external_refs_honours_an_xr_prefix(tmp_path: Path):
    """`\\externaldocument[M-]{main}` imports every label as `M-<label>`.

    The supplementary uses the prefix because xr also imports the paper's
    \\bibcite entries, which collide with the supplementary's own reference
    list; the references it makes to the paper are then `\\ref{M-fig:forest}`.
    """
    _external(tmp_path)
    src = "\\externaldocument[M-]{main}\nas Fig.~\\ref{M-fig:forest}.\n"
    assert "as Fig.~3." in bmd.resolve_external_refs(src, tmp_path)


def test_resolve_external_refs_refuses_an_unknown_label(tmp_path: Path):
    _external(tmp_path)
    src = "\\externaldocument{main}\nsee \\ref{fig:nowhere}\n"
    with pytest.raises(KeyError, match="fig:nowhere"):
        bmd.resolve_external_refs(src, tmp_path)


# -- figure captions above their graphics ------------------------------------
# Every figure caption in the manuscript sits above its graphic, which is what
# the journal asks for. Pandoc attaches a caption only when \caption follows
# the image, so read literally it drops all five without a word -- which is
# exactly what the first build after the move did.

FIGURE = (
    "\\begin{figure}[!ht]\n"
    "  \\caption{A caption with a brace group \\textbf{like this} in it.}\n"
    "  \\label{fig:one}\n"
    "  \\centering\n"
    "  \\includegraphics[width=\\linewidth]{figures/one.pdf}\n"
    "\\end{figure}\n"
)


def test_hoist_lifts_the_caption_in_front_of_the_float():
    out = bmd.hoist_figure_captions(FIGURE)
    assert out.index("A caption with a brace group") < out.index("\\begin{figure}")
    assert "\\caption{" not in out
    # The label stays with the float; only the caption text is lifted.
    assert "\\label{fig:one}" in out
    assert "\\includegraphics[width=\\linewidth]{figures/one.pdf}" in out


def test_hoist_keeps_a_nested_brace_group_whole():
    out = bmd.hoist_figure_captions(FIGURE)
    assert "\\textbf{like this} in it." in out
    assert out.count("\\begin{figure}") == 1
    assert out.count("\\end{figure}") == 1


def test_hoist_leaves_a_caption_that_already_follows_its_graphic():
    src = (
        "\\begin{figure}\n"
        "  \\includegraphics{figures/one.pdf}\n"
        "  \\caption{Pandoc reads this arrangement correctly.}\n"
        "\\end{figure}\n"
    )
    assert bmd.hoist_figure_captions(src) == src


def test_hoist_ignores_a_figure_named_inside_a_comment():
    # The preamble carries "\begin{figure*} for the volcano comparison" in a
    # comment explaining why that float is full-width. Without a line-start
    # anchor the scan begins there and swallows four figures into one match.
    src = "% so \\begin{figure*} for the volcano comparison\n" + FIGURE + FIGURE
    out = bmd.hoist_figure_captions(src)
    assert out.count("A caption with a brace group") == 2
    assert out.count("\\begin{figure}") == 2


def test_hoist_handles_the_starred_environment():
    src = FIGURE.replace("{figure}", "{figure*}")
    out = bmd.hoist_figure_captions(src)
    assert out.index("A caption with a brace group") < out.index("\\begin{figure*}")
    assert "\\end{figure*}" in out


# -- the abstract ---------------------------------------------------------------
# Pandoc lifts an abstract environment into document metadata and prints it
# straight after the author line, ahead of the affiliations, running title and
# correspondence the journal wants on the title page. As a headed section it
# stays where the source puts it.


def test_abstract_becomes_a_headed_section_in_place():
    src = "title page\n\\begin{abstract}\nThe text.\n\\end{abstract}\nKeywords\n"
    out = bmd.unwrap_abstract(src)
    assert "\\begin{abstract}" not in out and "\\end{abstract}" not in out
    assert out.index("title page") < out.index("\\section*{Abstract}") \
        < out.index("The text.") < out.index("Keywords")


def test_abstract_inside_twocolumn_is_still_recovered():
    """The Journal of Pain layout wrapped the abstract in \\twocolumn[...]."""
    src = ("\\twocolumn[\n  \\maketitle\n  \\begin{onecolabstract}\n"
           "  The text.\n  \\end{onecolabstract}\n]\nBody\n")
    out = bmd.unwrap_abstract(src)
    assert "\\section*{Abstract}" in out and "The text." in out
    assert "\\twocolumn" not in out


# -- AMA .bbl blocks and repository paths ---------------------------------------


def test_bbl_keeps_a_title_that_opens_with_a_protected_word():
    """Pandoc reads \\newblock as taking an argument and swallows the group.

    "{NCBI GEO}: archive for ..." reached the docx as ": archive for ...", and
    so did every other title opening with a braced acronym ({METAL}, {BDNF}).
    """
    bbl = ("\\begin{thebibliography}{1}\n\n\\bibitem{alpha}\nBarrett T.\n"
           "\\newblock {NCBI GEO}: archive for functional genomics data sets.\n"
           "\\newblock \\emph{Nucleic Acids Res}. 2013;41:D991-D995.\n\n"
           "\\end{thebibliography}\n")
    out = bmd.bbl_as_numbered_paragraphs(bbl, {"alpha": 1})
    assert "\\newblock" not in out
    assert "{NCBI GEO}: archive" in out


def test_repopath_reaches_the_docx_as_code():
    """Pandoc drops \\path, and \\repopath with it, leaving "reason codes ()"."""
    out = bmd.expand_repository_paths(
        "codes (\\repopath{literature/prisma/}) and \\repopath{scripts/a_b.py}.")
    assert out == ("codes (\\texttt{literature/prisma/}) and "
                   "\\texttt{scripts/a\\_b.py}.")


def test_native_numbering_follows_the_source():
    """Number captions in Word only where the source numbers them.

    The supplementary labels its one table "Supplementary Table 18" itself;
    pandoc's own numbering called it "Table 1".
    """
    assert bmd.numbers_its_captions("\\caption{A table.}")
    assert not bmd.numbers_its_captions("\\caption*{\\textbf{Supplementary Table 18} A.}")
