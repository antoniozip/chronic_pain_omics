"""Metabolomics Workbench discovery — the second metabolomics repository.

Until 2026-08-27 the metabolomic arm searched MetaboLights only. Every one of
the 39 records in `data/interim/metabolomics/harmonized.jsonl` carries
`source: METABOLIGHTS`, and `plan/cohort_expansion_plan.md:26` had already
flagged Metabolomics Workbench as unsearched. Three studies survived triage,
one of them carrying 90% of the pooled weight, so a second repository is the
only lever that changes the arm's shape rather than its numbers.

**Discovery enumerates the repository rather than querying it.** Workbench
exposes `/rest/study/study_title/<kw>/summary`, and it is not trustworthy for
systematic retrieval: `pain` returns a single study while `chronic fatigue`
returns fifteen, and truncated stems match inconsistently. `study_id/ST`
instead returns **all 4,507 studies** in one response (~85 s), so the whole
repository is retrieved and screened locally against one explicit expression.
That is the PRISMA order — retrieve with high sensitivity, exclude at triage
with a reason code — and it removes the query endpoint's behaviour from the
inclusion criteria entirely.

**Two response shapes, and the difference is a silent miscount.** A query
matching several studies returns an object keyed by rank, `{"1": {...},
"2": {...}}`. A query matching exactly one returns the **study object itself**,
whose ~15 metadata fields would be read as 15 studies by anything that counts
`values()`. This is the MetaboLights `/ws/studies` bare-string defect wearing a
different hat, and it is why `_studies()` tests for a `study_id` key before
unwrapping.

**Screening is on titles, because summaries carry no abstract.** The endpoint
returns title, species, institute, analysis type, sample count and dates —
enough to screen species and phenotype, not enough to judge a contrast. Studies
that pass here still face file triage, exactly as the MetaboLights nine did.
"""

from __future__ import annotations

import logging
import re
from typing import Any

from .base import DatasetRecord

logger = logging.getLogger(__name__)

SUMMARY_URL = "https://www.metabolomicsworkbench.org/rest/study/study_id/{prefix}/summary"
STUDY_UI = "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={acc}"

# Every accession in the repository shares this prefix, so it is also the
# enumerate-everything query.
STUDY_PREFIX = "ST"

# Pain vocabulary, as stems. Each alternative is anchored at its start only:
# a trailing \b after a truncated stem never matches, so `\barthrit\b` silently
# excludes "arthritis" — which cut a first screen of this repository from 23
# candidates to 7 before it was caught.
PAIN_PATTERN = re.compile(
    r"(?i)\b("
    r"pain|nocicept|noxious|analges|allodyn|hyperalges|neuralgi|migrain|"
    r"fibromyalg|arthrit|neuropath|sciatic|nerve injur|nerve ligat|"
    r"complete freund|endometrios|cystitis|irritable bowel|dysmenorrh|"
    r"headache|crps|causalgi|chemotherapy.induced periph|cipn|low back|"
    r"disc degener|trigemin|temporomandib"
    r")"
)

# Species accepted by the study's inclusion criteria. Workbench reports one
# species string per study, unlike MetaboLights' list.
SPECIES_KEEP = {"homo sapiens", "mus musculus", "rattus norvegicus"}

# "lipidomic" (no s) is a prefix of "lipidomics", matching the single entry
# MetabolomicsIngester uses for the same purpose.
LIPID_KEYWORDS = ("lipidomic", "lipidome", "lipid profiling")


def _studies(payload: Any) -> list[dict]:
    """Normalise Workbench's two response shapes to a list of study objects.

    A single-hit response is the bare study object; counting its fields would
    report one study as fifteen. See the module docstring.
    """
    if not payload:
        return []
    if isinstance(payload, dict):
        if "study_id" in payload:
            return [payload]
        return [v for v in payload.values() if isinstance(v, dict) and "study_id" in v]
    if isinstance(payload, list):
        return [v for v in payload if isinstance(v, dict) and "study_id" in v]
    return []


def is_pain_related(title: str) -> bool:
    """Whether a study title carries pain vocabulary."""
    return bool(PAIN_PATTERN.search(title or ""))


def is_lipidomics(text: str) -> bool:
    return any(kw in (text or "").lower() for kw in LIPID_KEYWORDS)


def species_accepted(species: str) -> bool:
    return (species or "").strip().lower() in SPECIES_KEEP


def _n_samples(raw: Any) -> int:
    """Sample count, or -1 for unknown.

    Workbench reports this as a string. -1 rather than 0 because
    `apply_filters` passes -1 through as unknown but reads 0 as "fewer than
    min_samples" and drops the study — the same convention the MetaboLights
    ingester uses for a field EBI Search does not carry at all.
    """
    try:
        n = int(str(raw).strip())
    except (TypeError, ValueError):
        return -1
    return n if n > 0 else -1


def to_record(study: dict, modality: str) -> DatasetRecord:
    """Build a DatasetRecord from one Workbench study summary.

    The URL is constructed from the accession rather than taken from the
    payload's `study_url`, which echoes the query keyword
    (`...StudyID=ARTHRITIS`) instead of naming the study.
    """
    acc = study.get("study_id", "")
    species = (study.get("species") or "").strip()
    return DatasetRecord(
        modality=modality,
        source="METABOLOMICS_WORKBENCH",
        accession=acc,
        title=(study.get("study_title") or "").strip(),
        # Summaries carry no abstract; institute and analysis type are what
        # there is, and screening therefore runs on the title.
        description=(study.get("institute") or "").strip(),
        species=[species] if species else [],
        n_samples=_n_samples(study.get("number_of_samples")),
        platform=(study.get("analysis_type") or "Metabolomics").strip(),
        year=(study.get("submission_date") or "")[:4],
        url=STUDY_UI.format(acc=acc),
        keywords=[],
    )


# Enumerating all 4,507 studies takes ~85 s. The base class's 30 s default
# would time out on every call, so this is passed explicitly rather than left
# to DEFAULT_TIMEOUT.
ENUMERATE_TIMEOUT = 300


def fetch_all(getter: Any, prefix: str = STUDY_PREFIX) -> list[dict]:
    """Enumerate the repository. `getter` is OmicsIngester._get."""
    payload = getter(SUMMARY_URL.format(prefix=prefix), timeout=ENUMERATE_TIMEOUT)
    studies = _studies(payload)
    logger.info("[workbench] enumerated %d studies", len(studies))
    return studies


def screen(studies: list[dict], modality: str = "metabolomics") -> tuple[list[dict], list[dict]]:
    """Split enumerated studies into kept and excluded, with a reason code each.

    Exclusions are returned rather than dropped so every one reaches the PRISMA
    record. Reason codes match those already in
    `literature/prisma/metabolomics_manual_review.csv`.
    """
    kept: list[dict] = []
    excluded: list[dict] = []
    want_lipids = modality == "lipidomics"
    for s in studies:
        acc = s.get("study_id", "")
        title = s.get("study_title", "") or ""
        species = s.get("species", "") or ""
        row = {"accession": acc, "title": title, "species": species,
               "n_samples": s.get("number_of_samples", "")}
        if not is_pain_related(title):
            excluded.append({**row, "reason_code": "no_pain_phenotype"})
            continue
        if not species_accepted(species):
            excluded.append({**row, "reason_code": "wrong_species"})
            continue
        # A study belongs to exactly one of the two arms, so the lipidomics
        # ingester inverts this test rather than keeping its own keyword list.
        if is_lipidomics(title) != want_lipids:
            excluded.append({**row, "reason_code": "wrong_modality_arm"})
            continue
        kept.append(s)
    return kept, excluded
