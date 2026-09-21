"""Single source of truth for pain-model stratum names and their labels.

06b writes one `{stem}_pooled.csv` per stratum, and the stems it chooses
(`other_rodent_amendment`, `neuropathy_rodent`) are file identifiers rather
than anything a reader should see. Two places turn them into prose: the
stratified table in the manuscript and panel B of the graphical abstract.

They used to hold separate lists, which meant the figure could call a stratum
`other_rodent_amendment` while the table beside it called the same stratum
"Other (rodent)" -- the state this module was created to end. A reader
comparing the two has no way to tell that they are the same row.

Labels here are plain text, since that is what matplotlib renders and what
LaTeX accepts unescaped for all but one of them. `LATEX_LABELS` carries the
exceptions, so markup lives with the writer that needs it rather than in the
data.
"""

from __future__ import annotations

#: (file stem, printed label) for the strata of the original design.
#: Order is the printed order.
ORIGINAL: list[tuple[str, str]] = [
    ("CCI", "CCI"),
    ("SNI", "SNI"),
    ("SNL", "SNL"),
    ("CFA", "CFA"),
    ("CIPN", "CIPN"),
    ("neuropathic_human", "Neuropathic (human)"),
    ("lbp_human", "LBP (human)"),
    ("nociplastic_human", "Nociplastic (human)"),
    ("invitro_human", "In vitro (human)"),
]

#: Strata added by the search amendment, reported below the rule in the table.
AMENDMENT: list[tuple[str, str]] = [
    ("endometriosis_human", "Endometriosis (human)"),
    ("endometriosis_rodent", "Endometriosis (rodent)"),
    ("fibromyalgia_human", "Fibromyalgia (human)"),
    ("ibs_human", "IBS (human)"),
    ("neuropathy_human", "Neuropathy (human)"),
    ("neuropathy_rodent", "Neuropathy (rodent)"),
    ("other_human_amendment", "Other (human)"),
    ("other_rodent_amendment", "Other (rodent)"),
]

#: Stem -> label, for callers that do not care about print order.
STRATUM_LABELS: dict[str, str] = dict(ORIGINAL + AMENDMENT)

#: Labels needing LaTeX markup the plain form cannot carry. Applied by the
#: table writer; every other stratum prints its plain label unchanged.
LATEX_LABELS: dict[str, str] = {
    "invitro_human": r"\emph{In vitro} (human)",
}


def latex_label(stem: str) -> str:
    """Return the stratum's label with LaTeX markup where it has any."""
    return LATEX_LABELS.get(stem, STRATUM_LABELS[stem])


def plain_label(stem: str) -> str:
    """Return the stratum's label as plain text, for figures and logs.

    Unknown stems return themselves: a figure should draw the stem rather
    than raise, the table being the place that insists on completeness.
    """
    return STRATUM_LABELS.get(stem, stem)


__all__ = [
    "ORIGINAL",
    "AMENDMENT",
    "STRATUM_LABELS",
    "LATEX_LABELS",
    "latex_label",
    "plain_label",
]
