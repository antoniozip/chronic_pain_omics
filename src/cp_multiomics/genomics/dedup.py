"""Cohort-independence dedup for GWAS summary-statistics studies.

Implements plan/cohort_expansion_plan.md §7.7: inverse-variance meta-analysis
assumes independent samples, but the retrieved pain studies overlap heavily
(UK Biobank re-analyses) and include gene-based burden analyses that are not
SNP-level GWAS. This reduces the ~138 retrievable studies to the independent,
poolable subset BEFORE any large download.
"""

from __future__ import annotations

import dataclasses
import logging

from cp_multiomics.genomics.resolve import SumstatsStudy

logger = logging.getLogger(__name__)

PRESUMED_BIOBANK = "PRESUMED_BIOBANK"
BIOBANK_N_THRESHOLD = 100_000


def apply_presumed_cohort(
    studies: list[SumstatsStudy],
    n_threshold: int = BIOBANK_N_THRESHOLD,
) -> list[SumstatsStudy]:
    """Reassign empty-cohort biobank-scale / ICD10-coded studies to a shared
    presumed-biobank cohort, so §7.7 collapses their duplicate phenotypes.

    detect_cohort is text-only and misses UK Biobank re-analyses whose sample
    text omits the literal cohort name; without this, large-N re-analyses of the
    same biobank are wrongly treated as independent, under-collapsing the cohort
    and inflating meta-analysis false significance. Errs toward independence.

    A study is reassigned iff its cohort is currently "" AND (n > n_threshold OR
    its phenotype is ICD10-coded). Studies with a real cohort, or small-N
    non-ICD10 empty-cohort studies, are returned unchanged.
    """
    out: list[SumstatsStudy] = []
    for s in studies:
        if s.cohort == "" and (s.n > n_threshold or s.phenotype.startswith("icd10")):
            out.append(dataclasses.replace(s, cohort=PRESUMED_BIOBANK))
        else:
            out.append(s)
    return out


# PRESUMED_BIOBANK studies are inferred UK Biobank (large-N or ICD10-coded, no
# literal cohort name). For dedup grouping they must share a bucket with studies
# whose cohort was literally detected as UKB, so same-phenotype re-analyses of UK
# Biobank collapse. The study keeps its original label in the manifest (so the
# named-vs-inferred split stays auditable); only the GROUP KEY is canonicalized.
_COHORT_ALIASES = {PRESUMED_BIOBANK: "UKB"}


def _canonical_cohort(cohort: str) -> str:
    return _COHORT_ALIASES.get(cohort, cohort)


def dedup_studies(
    studies: list[SumstatsStudy],
) -> tuple[list[SumstatsStudy], list[tuple[SumstatsStudy, str]]]:
    """Return (kept, dropped) applying the §7.7 independence rule.

    - Gene-based burden studies are dropped (reason "burden").
    - Among studies sharing a non-empty (cohort, phenotype), keep the largest N
      (ties broken by lexicographically smallest accession); the rest are
      dropped (reason "cohort_overlap").
    - Studies with an empty cohort are presumed independent and never collapsed.
    """
    kept: list[SumstatsStudy] = []
    dropped: list[tuple[SumstatsStudy, str]] = []

    non_burden: list[SumstatsStudy] = []
    for s in studies:
        if s.is_burden:
            dropped.append((s, "burden"))
        else:
            non_burden.append(s)

    groups: dict[tuple[str, str], list[SumstatsStudy]] = {}
    for s in non_burden:
        if s.cohort == "":
            kept.append(s)                       # presumed independent
            continue
        groups.setdefault((_canonical_cohort(s.cohort), s.phenotype), []).append(s)

    for members in groups.values():
        max_n = max(s.n for s in members)
        candidates = [s for s in members if s.n == max_n]
        winner = min(candidates, key=lambda s: s.accession)
        kept.append(winner)
        for s in members:
            if s.accession != winner.accession:
                dropped.append((s, "cohort_overlap"))

    logger.info("Dedup: %d kept, %d dropped (%d burden, %d overlap)",
                len(kept), len(dropped),
                sum(1 for _, r in dropped if r == "burden"),
                sum(1 for _, r in dropped if r == "cohort_overlap"))
    return kept, dropped
