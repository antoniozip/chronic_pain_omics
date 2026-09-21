"""Resolve GWAS Catalog studies that carry full summary statistics.

The per-study summary-statistics REST endpoint 404s and the HAL
`summaryStatistics` link is empty (verified 2026-07-11), so studies are
discovered through the trait-search REST API filtered on `fullPvalueSet`, and
the sumstats files themselves are located via the FTP path convention (see
cp_multiomics.genomics.ftp).
"""

from __future__ import annotations

import logging
import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from cp_multiomics.genomics.sample_size import (
    detect_cohort,
    is_burden_study,
    parse_sample_size,
)

logger = logging.getLogger(__name__)

GWAS_REST = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByEfoTrait"

_PAREN_RE = re.compile(r"\(.*?\)")
_WS_RE = re.compile(r"\s+")


@dataclass(frozen=True)
class SumstatsStudy:
    accession: str
    trait: str
    sample_size_text: str
    n: int
    cohort: str
    is_burden: bool
    phenotype: str


def normalize_phenotype(trait: str) -> str:
    """Lowercased, parenthetical-stripped, whitespace-collapsed trait key.

    If stripping parentheticals empties the string (a trait that is ENTIRELY a
    parenthetical annotation, e.g. "(UKB data field 3773)"), fall back to the
    normalized full trait rather than "" — an empty key would collide unrelated
    studies in the (cohort, phenotype) dedup and silently drop real studies.
    """
    stripped = _WS_RE.sub(" ", _PAREN_RE.sub("", trait or "")).strip().lower()
    if stripped:
        return stripped
    # Trait is entirely a parenthetical annotation (e.g. "(UKB data field
    # 3773)"): strip only the paren characters, keeping the inner content,
    # instead of collapsing to "".
    unwrapped = re.sub(r"[()]", "", trait or "")
    return _WS_RE.sub(" ", unwrapped).strip().lower()


def resolve_sumstats_studies(
    traits: list[str],
    http_get: Callable[..., dict[str, Any]],
) -> list[SumstatsStudy]:
    """Return one SumstatsStudy per unique accession with fullPvalueSet=True.

    Args:
        traits: EFO trait strings to query.
        http_get: Callable(url, params=...) -> parsed JSON dict. Injected so the
            suite can run without network.

    Returns:
        Deduplicated list of SumstatsStudy (first occurrence wins).
    """
    seen: dict[str, SumstatsStudy] = {}
    for trait in traits:
        data = http_get(GWAS_REST, params={"efoTrait": trait, "size": 500})
        studies = data.get("_embedded", {}).get("studies", []) or []
        for raw in studies:
            if not raw.get("fullPvalueSet"):
                continue
            acc = raw.get("accessionId", "")
            if not acc or acc in seen:
                continue
            trait_str = (raw.get("diseaseTrait") or {}).get("trait", "") or ""
            sample_text = raw.get("initialSampleSize", "") or ""
            seen[acc] = SumstatsStudy(
                accession=acc,
                trait=trait_str,
                sample_size_text=sample_text,
                n=parse_sample_size(sample_text),
                cohort=detect_cohort(f"{sample_text} {trait_str}"),
                is_burden=is_burden_study(trait_str),
                phenotype=normalize_phenotype(trait_str),
            )
    logger.info("Resolved %d unique sumstats studies", len(seen))
    return list(seen.values())
