"""Omics dataset ingestion package.

Each modality is registered via @register_ingester('name') and
instantiated through IngesterFactory(modality).
"""

from __future__ import annotations

from .base import DatasetRecord, OmicsIngester

INGESTER_REGISTRY: dict[str, type[OmicsIngester]] = {}


def register_ingester(name: str):
    """Class decorator that registers an OmicsIngester under a modality name."""
    def decorator(cls: type[OmicsIngester]) -> type[OmicsIngester]:
        INGESTER_REGISTRY[name] = cls
        return cls
    return decorator


def IngesterFactory(modality: str) -> OmicsIngester:
    """Instantiate the registered ingester for a given modality."""
    if modality not in INGESTER_REGISTRY:
        available = sorted(INGESTER_REGISTRY)
        raise ValueError(f"Unknown modality '{modality}'. Available: {available}")
    return INGESTER_REGISTRY[modality]()


# Import all ingesters to trigger registration
from .genomics import GenomicsIngester  # noqa: E402, F401
from .lipidomics import LipidomicsIngester  # noqa: E402, F401
from .metabolomics import MetabolomicsIngester  # noqa: E402, F401
from .proteomics import ProteomicsIngester  # noqa: E402, F401
from .transcriptomics import TranscriptomicsIngester  # noqa: E402, F401

__all__ = [
    "DatasetRecord",
    "OmicsIngester",
    "IngesterFactory",
    "register_ingester",
    "INGESTER_REGISTRY",
]
