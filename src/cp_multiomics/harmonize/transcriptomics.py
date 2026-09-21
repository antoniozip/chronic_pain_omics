"""Transcriptomics harmonizer — GEO datasets (RNA-seq, microarray)."""

from __future__ import annotations

import logging
from typing import Any

from . import register_harmonizer
from .base import HarmonizedRecord, OmicsHarmonizer

logger = logging.getLogger(__name__)

# GEO platform accessions for common RNA-seq platforms
RNASEQ_GPL_PREFIXES = {"GPL24676", "GPL20301", "GPL18573", "GPL16791", "GPL21290"}


@register_harmonizer("transcriptomics")
class TranscriptomicsHarmonizer(OmicsHarmonizer):
    """Harmonizes GEO transcriptomics metadata to canonical platform/species/ID space."""

    modality = "transcriptomics"

    def _harmonize_one(self, raw: dict[str, Any], cfg: dict[str, Any]) -> HarmonizedRecord:
        min_samples = cfg.get("filters", {}).get("min_samples", 5)

        species_raw = raw.get("species", [])
        species_canonical = self.normalize_species(species_raw)
        has_human = "Homo sapiens" in species_canonical
        has_animal = any(s in species_canonical for s in ("Mus musculus", "Rattus norvegicus"))

        platform_raw = raw.get("platform", "")
        platform_canonical = self.normalize_platform(platform_raw)

        # Determine expected ID space from platform
        id_space = self._infer_id_space(platform_canonical, raw.get("keywords", []))

        # For microarray, effect size is typically log2FC; for RNA-seq also log2FC or VST
        effect_size_unit = "log2FC"

        return HarmonizedRecord(
            modality=self.modality,
            source=raw.get("source", "GEO"),
            accession=raw.get("accession", ""),
            title=raw.get("title", ""),
            species_canonical=species_canonical,
            has_human=has_human,
            has_animal=has_animal,
            platform_canonical=platform_canonical,
            year=raw.get("year", ""),
            n_samples=int(raw.get("n_samples", 0)),
            id_space=id_space,
            effect_size_unit=effect_size_unit,
            url=raw.get("url", ""),
            ortholog_needed=has_animal and not has_human,
            quality_flags=self.flag_quality(raw, min_samples),
        )

    @staticmethod
    def _infer_id_space(platform: str, keywords: list[str]) -> str:
        kw_text = " ".join(keywords).lower()
        if "ensembl" in kw_text:
            return "ENSEMBL"
        if "entrez" in kw_text:
            return "ENTREZ"
        # Default for GEO: HGNC symbols after remapping
        return "HGNC"
