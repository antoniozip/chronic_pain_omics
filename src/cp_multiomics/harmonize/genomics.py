"""Genomics harmonizer — GWAS Catalog studies."""

from __future__ import annotations

import logging
from typing import Any

from . import register_harmonizer
from .base import HarmonizedRecord, OmicsHarmonizer

logger = logging.getLogger(__name__)


@register_harmonizer("genomics")
class GenomicsHarmonizer(OmicsHarmonizer):
    """Harmonizes GWAS Catalog metadata.

    GWAS Catalog studies are always human-only.
    Effect size unit defaults to 'beta' (continuous traits) or 'OR' (binary).
    ID space is rsID for variants.
    Ortholog mapping is not applicable.
    """

    modality = "genomics"

    def _harmonize_one(self, raw: dict[str, Any], cfg: dict[str, Any]) -> HarmonizedRecord:
        min_samples = cfg.get("filters", {}).get("min_samples", 5)

        # GWAS Catalog is human-only
        species_canonical = ["Homo sapiens"]

        platform_canonical = "GWAS"

        # Infer effect-size unit from study description / keywords
        desc = (raw.get("description", "") + " " + " ".join(raw.get("keywords", []))).lower()
        effect_size_unit = self._infer_effect_size(desc)

        flags = self.flag_quality(raw, min_samples)
        if int(raw.get("n_samples", 0)) < 1000:
            flags.append("small_gwas(<1000)")

        return HarmonizedRecord(
            modality=self.modality,
            source=raw.get("source", "GWAS_CATALOG"),
            accession=raw.get("accession", ""),
            title=raw.get("title", ""),
            species_canonical=species_canonical,
            has_human=True,
            has_animal=False,
            platform_canonical=platform_canonical,
            year=raw.get("year", ""),
            n_samples=int(raw.get("n_samples", 0)),
            id_space="rsID",
            effect_size_unit=effect_size_unit,
            url=raw.get("url", ""),
            ortholog_needed=False,
            quality_flags=flags,
        )

    @staticmethod
    def _infer_effect_size(description_lower: str) -> str:
        if "odds ratio" in description_lower or " or " in description_lower:
            return "OR"
        if "hazard ratio" in description_lower or " hr " in description_lower:
            return "HR"
        # Default for quantitative pain traits (NRS, VAS, etc.)
        return "beta"
