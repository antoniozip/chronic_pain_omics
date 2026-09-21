"""Omics harmonization package.

Each modality has a registered Harmonizer that standardizes identifiers,
units, and metadata from the raw ingest manifests (data/raw/{modality}/datasets.jsonl)
into harmonized tables (data/interim/{modality}/harmonized.jsonl).

Usage pattern:
    harmonizer = HarmonizerFactory("transcriptomics")
    records = harmonizer.harmonize(raw_records, cfg)
"""

from __future__ import annotations

from .base import HarmonizedRecord, OmicsHarmonizer

HARMONIZER_REGISTRY: dict[str, type[OmicsHarmonizer]] = {}


def register_harmonizer(name: str):
    """Class decorator registering a harmonizer under a modality name."""
    def decorator(cls: type[OmicsHarmonizer]) -> type[OmicsHarmonizer]:
        HARMONIZER_REGISTRY[name] = cls
        return cls
    return decorator


def HarmonizerFactory(modality: str) -> OmicsHarmonizer:
    if modality not in HARMONIZER_REGISTRY:
        available = sorted(HARMONIZER_REGISTRY)
        raise ValueError(f"Unknown modality '{modality}'. Available: {available}")
    return HARMONIZER_REGISTRY[modality]()


from .genomics import GenomicsHarmonizer  # noqa: E402, F401
from .lipidomics import LipidomicsHarmonizer  # noqa: E402, F401
from .metabolomics import MetabolomicsHarmonizer  # noqa: E402, F401
from .proteomics import ProteomicsHarmonizer  # noqa: E402, F401
from .transcriptomics import TranscriptomicsHarmonizer  # noqa: E402, F401

__all__ = [
    "HarmonizedRecord",
    "OmicsHarmonizer",
    "HarmonizerFactory",
    "register_harmonizer",
    "HARMONIZER_REGISTRY",
]
