"""Parse Metabolomics Workbench `mwtab` files into abundance matrices.

Workbench deposits one text file per study holding metadata, the abundance
matrix and the metabolite annotation. It is not MetaboLights' ISA-Tab, and the
differences that matter are these.

**One file can hold several analyses.** A study measured in positive and
negative ion mode concatenates two complete documents, each opening with
`#METABOLOMICS WORKBENCH ... ANALYSIS_ID:ANxxxxxx`. These are the *same
samples* measured twice, exactly like MetaboLights' multiple MAFs, so they are
parsed separately and collapsed per study afterwards. Treating them as
independent units is the pseudo-replication that let PXD013362 enter the
proteomic pool fourteen times.

**Sample factors are inline.** The data block's second row carries each
column's full factor string, so mapping a column to a study arm needs no join
against a sample table — none of the three-hop resolution MetaboLights
requires.

**Below-detection values are not zero.** Workbench writes them as `< LOD`,
`< LLOQ`, `NA`, or an empty cell, and the comparison operators arrive
HTML-escaped as `&lt;`. ST003984 has 6,300 `&lt; LOD` cells and ST003177 has
41,245 `NA`. Read as zero, a metabolite absent from every control and present
in every case yields an enormous fabricated effect size; read as missing, it
simply carries less weight. `parse_value` therefore returns NaN for all of
them, and never 0.0.
"""

from __future__ import annotations

import html
import logging
import re
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

# Opens each analysis document within an mwtab file.
ANALYSIS_HEADER = re.compile(r"^#METABOLOMICS WORKBENCH .*?ANALYSIS_ID:(AN\d+)", re.M)

# The abundance matrix. NMR studies use a differently-named block, so both are
# recognised; a study deposited only as binned NMR carries no named
# metabolites and is excluded upstream at triage rather than here.
DATA_START = re.compile(r"^(MS_METABOLITE_DATA|NMR_METABOLITE_DATA|NMR_BINNED_DATA)_START",
                        re.M)
DATA_END = re.compile(r"^(MS_METABOLITE_DATA|NMR_METABOLITE_DATA|NMR_BINNED_DATA)_END",
                      re.M)
METABOLITES_START = re.compile(r"^METABOLITES_START", re.M)
METABOLITES_END = re.compile(r"^METABOLITES_END", re.M)
UNITS = re.compile(r"^\w+_METABOLITE_DATA:UNITS\t(.+)$", re.M)

# Cells that mean "not measurable here". Matched after HTML-unescaping, so
# `&lt; LOD` has already become `< LOD`.
MISSING = re.compile(r"^\s*(|na|n/?a|nd|nan|null|-|\.|<\s*\w+|below\s.*)\s*$", re.I)


@dataclass
class Analysis:
    """One analysis document: its samples, their factors, and the matrix."""

    analysis_id: str
    units: str = ""
    samples: list[str] = field(default_factory=list)
    factors: list[str] = field(default_factory=list)
    # metabolite name -> one value per sample, aligned with `samples`
    values: dict[str, list[float]] = field(default_factory=dict)
    # metabolite name -> RefMet name, where the study supplies one
    refmet: dict[str, str] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(self.samples)


def parse_value(raw: str) -> float:
    """One abundance cell as a float, or NaN where it is not a measurement.

    Never returns 0.0 for a below-detection marker. `< LOD` means the compound
    was not measurable, not that its abundance was zero, and substituting zero
    manufactures effect sizes in exactly the comparisons the arm exists to
    make.
    """
    text = html.unescape(raw or "").strip()
    if MISSING.match(text):
        return float("nan")
    try:
        return float(text.replace(",", ""))
    except ValueError:
        return float("nan")


def split_analyses(text: str) -> list[tuple[str, str]]:
    """Split an mwtab file into (analysis_id, document) pairs.

    A file with one analysis yields one pair. A file with none yields none,
    rather than treating the whole text as a single untitled analysis, so a
    truncated download cannot be mistaken for data.
    """
    matches = list(ANALYSIS_HEADER.finditer(text))
    if not matches:
        return []
    out = []
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        out.append((m.group(1), text[m.start():end]))
    return out


def _block(text: str, start: re.Pattern, end: re.Pattern) -> list[str]:
    s = start.search(text)
    if not s:
        return []
    e = end.search(text, s.end())
    stop = e.start() if e else len(text)
    return [ln for ln in text[s.end():stop].split("\n") if ln.strip()]


def refmet_map(document: str) -> dict[str, str]:
    """`metabolite_name -> RefMet Name` from the METABOLITES block.

    RefMet is a standardised naming space, so it resolves to ChEBI far more
    reliably than the free-text name a study happens to use. Absent for older
    studies, in which case the raw name is all there is.
    """
    lines = _block(document, METABOLITES_START, METABOLITES_END)
    if not lines:
        return {}
    header = [c.strip().lower() for c in lines[0].split("\t")]
    try:
        name_i = header.index("metabolite_name")
    except ValueError:
        return {}
    refmet_i = next((i for i, c in enumerate(header) if c.startswith("refmet")), -1)
    if refmet_i < 0:
        return {}
    out = {}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) > max(name_i, refmet_i):
            name, refmet = cells[name_i].strip(), cells[refmet_i].strip()
            # "-" is Workbench's "no RefMet assignment", not a name.
            if name and refmet and refmet != "-":
                out[name] = refmet
    return out


def parse_analysis(analysis_id: str, document: str) -> Analysis | None:
    """One analysis document into an Analysis, or None if it holds no matrix."""
    lines = _block(document, DATA_START, DATA_END)
    if len(lines) < 3:
        return None

    header = lines[0].split("\t")
    if header[0].strip().lower() != "samples":
        logger.warning("%s: data block does not open with a Samples row", analysis_id)
        return None
    samples = [c.strip() for c in header[1:]]

    factors: list[str] = []
    first_metabolite = 1
    if lines[1].split("\t")[0].strip().lower() == "factors":
        factors = [c.strip() for c in lines[1].split("\t")[1:]]
        first_metabolite = 2

    # A short Factors row would silently misalign every column against the
    # wrong arm, so the two are padded to a common length rather than zipped.
    if factors and len(factors) != len(samples):
        logger.warning("%s: %d samples but %d factor cells; padding",
                       analysis_id, len(samples), len(factors))
        factors = (factors + [""] * len(samples))[:len(samples)]

    units_match = UNITS.search(document)
    out = Analysis(analysis_id=analysis_id,
                   units=units_match.group(1).strip() if units_match else "",
                   samples=samples, factors=factors,
                   refmet=refmet_map(document))

    for line in lines[first_metabolite:]:
        cells = line.split("\t")
        name = cells[0].strip()
        if not name:
            continue
        row = [parse_value(c) for c in cells[1:]]
        # Pad or trim to the sample count for the same alignment reason.
        row = (row + [float("nan")] * len(samples))[:len(samples)]
        if np.isfinite(row).any():
            out.values[name] = row
    return out


def parse_mwtab(text: str) -> list[Analysis]:
    """Every analysis in an mwtab file that carries an abundance matrix."""
    out = []
    for analysis_id, document in split_analyses(text):
        parsed = parse_analysis(analysis_id, document)
        if parsed is None:
            logger.info("  %s: no abundance matrix, skipped", analysis_id)
            continue
        out.append(parsed)
    return out
