"""Lipidomics harmonizer — MetaboLights (lipidomics-flagged) datasets."""

from __future__ import annotations

import logging
from typing import Any

from . import register_harmonizer
from .base import HarmonizedRecord, OmicsHarmonizer

logger = logging.getLogger(__name__)


@register_harmonizer("lipidomics")
class LipidomicsHarmonizer(OmicsHarmonizer):
    """Harmonizes MetaboLights lipidomics metadata.

    ID space: LIPID MAPS (LMID) as canonical; HMDB as fallback.
    Effect size unit: SMD across lipid concentrations.
    Platform is almost always LC-MS for modern lipidomics.
    """

    modality = "lipidomics"

    def _harmonize_one(self, raw: dict[str, Any], cfg: dict[str, Any]) -> HarmonizedRecord:
        min_samples = cfg.get("filters", {}).get("min_samples", 5)

        species_raw = raw.get("species", [])
        species_canonical = self.normalize_species(species_raw)
        has_human = "Homo sapiens" in species_canonical
        has_animal = any(s in species_canonical for s in ("Mus musculus", "Rattus norvegicus"))

        platform_raw = raw.get("platform", "")
        platform_canonical = self.normalize_platform(platform_raw)
        if platform_canonical == "unknown":
            platform_canonical = "LC-MS (lipidomics)"

        flags = self.flag_quality(raw, min_samples)

        return HarmonizedRecord(
            modality=self.modality,
            source=raw.get("source", "METABOLIGHTS"),
            accession=raw.get("accession", ""),
            title=raw.get("title", ""),
            species_canonical=species_canonical,
            has_human=has_human,
            has_animal=has_animal,
            platform_canonical=platform_canonical,
            year=raw.get("year", ""),
            n_samples=int(raw.get("n_samples", 0)),
            id_space="LIPIDMAPS",
            effect_size_unit="SMD",
            url=raw.get("url", ""),
            ortholog_needed=has_animal and not has_human,
            quality_flags=flags,
        )
