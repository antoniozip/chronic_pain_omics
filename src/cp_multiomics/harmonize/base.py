"""Base class and shared data model for all omics harmonizers.

Harmonization scope for a meta-study (metadata-level, not full raw data):
  - Canonical identifier assignment (HGNC symbol, UniProt, ChEBI, LIPID MAPS)
  - Technology / platform normalization to a controlled vocabulary
  - Species tag normalization to NCBI Taxonomy names
  - Ortholog flag for cross-species studies
  - Effect-size unit annotation (what unit will step 05 expect)
  - Study-level quality flags

Harmonizers do NOT download raw count matrices or full datasets —
that belongs in the per-study analysis layer (step 05).
"""

from __future__ import annotations

import json
import logging
import re
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# Controlled vocabulary for technology
PLATFORM_VOCAB: dict[str, str] = {
    # Transcriptomics
    "rna-seq": "RNA-seq",
    "rnaseq": "RNA-seq",
    "rna seq": "RNA-seq",
    "high throughput sequencing": "RNA-seq",
    "expression profiling by high throughput sequencing": "RNA-seq",
    "microarray": "microarray",
    "expression profiling by array": "microarray",
    "affymetrix": "microarray",
    "illumina": "RNA-seq",
    # Proteomics
    "mass spectrometry": "mass spectrometry",
    "lcms": "LC-MS",
    "lc-ms": "LC-MS",
    "lc ms": "LC-MS",
    # Metabolomics / Lipidomics
    "nmr": "NMR",
    "nmr spectroscopy": "NMR",
    "gcms": "GC-MS",
    "gc-ms": "GC-MS",
    "gc ms": "GC-MS",
    "lipidomics ms": "LC-MS (lipidomics)",
    # Genomics
    "gwas": "GWAS",
    "genome-wide association": "GWAS",
}

SPECIES_VOCAB: dict[str, str] = {
    "homo sapiens": "Homo sapiens",
    "human": "Homo sapiens",
    "mus musculus": "Mus musculus",
    "mouse": "Mus musculus",
    "rattus norvegicus": "Rattus norvegicus",
    "rat": "Rattus norvegicus",
}


@dataclass
class HarmonizedRecord:
    """Metadata-level harmonized record for one dataset."""
    modality: str
    source: str
    accession: str
    title: str
    species_canonical: list[str]        # NCBI Taxonomy names
    has_human: bool
    has_animal: bool
    platform_canonical: str             # controlled-vocabulary platform
    year: str
    n_samples: int
    # "HGNC" | "UniProt" | "ChEBI" | "LIPIDMAPS" | "rsID" | "unknown"
    id_space: str
    effect_size_unit: str               # unit expected in step 05 (e.g. "log2FC", "beta", "SMD")
    url: str
    ortholog_needed: bool               # True if non-human species present
    quality_flags: list[str] = field(default_factory=list)
    harmonized_at: str = field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self))


class OmicsHarmonizer(ABC):
    """Abstract base for per-modality harmonizers."""

    modality: str = ""

    def harmonize(
        self, raw_records: list[dict[str, Any]], cfg: dict[str, Any]
    ) -> list[HarmonizedRecord]:
        """Harmonize a list of raw DatasetRecord dicts into HarmonizedRecords."""
        results: list[HarmonizedRecord] = []
        failures: list[str] = []
        for raw in raw_records:
            try:
                rec = self._harmonize_one(raw, cfg)
                results.append(rec)
            except Exception as exc:
                failures.append(
                    f"{raw.get('source', '?')}/{raw.get('accession', '?')}: {exc}"
                )
                logger.warning(
                    "Failed to harmonize %s/%s: %s",
                    raw.get("source", "?"), raw.get("accession", "?"), exc,
                )
        if failures:
            logger.warning(
                "[%s] %d / %d records failed harmonization; first failures: %s",
                self.modality, len(failures), len(raw_records),
                "; ".join(failures[:5]),
            )
        logger.info("[%s] Harmonized %d / %d records",
                    self.modality, len(results), len(raw_records))
        if raw_records and not results:
            # An empty harmonized.jsonl is read downstream as "0 studies", the
            # same shape as a genuinely empty cohort. A modality that produced
            # records and harmonized none is a failure, not a result.
            raise RuntimeError(
                f"[{self.modality}] harmonized 0 of {len(raw_records)} records; "
                "refusing to write an empty manifest. First failures: "
                + "; ".join(failures[:5])
            )
        return results

    @abstractmethod
    def _harmonize_one(self, raw: dict[str, Any], cfg: dict[str, Any]) -> HarmonizedRecord:
        """Harmonize a single raw record dict."""

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------

    @staticmethod
    def normalize_platform(platform_raw: str) -> str:
        key = platform_raw.lower().strip()
        for pattern, canonical in PLATFORM_VOCAB.items():
            if pattern in key:
                return canonical
        return platform_raw or "unknown"

    @staticmethod
    def normalize_species(species_list: list[str]) -> list[str]:
        normalized: list[str] = []
        for s in species_list:
            # Strip parenthetical annotations e.g. "Mus musculus (mouse)" → "Mus musculus"
            key = re.sub(r"\s*\(.*?\)", "", s).lower().strip()
            normalized.append(SPECIES_VOCAB.get(key, s))
        return list(dict.fromkeys(normalized))  # deduplicate, preserve order

    @staticmethod
    def flag_quality(raw: dict[str, Any], min_samples: int) -> list[str]:
        flags: list[str] = []
        if raw.get("n_samples", 0) < min_samples:
            flags.append(f"low_n_samples(<{min_samples})")
        # Coerce before slicing: a provider that returns an integer year used
        # to raise TypeError here, and the broad except in harmonize() then
        # dropped the whole record without a trace.
        year = str(raw.get("year") or "")
        if not year or year[:4] == "0000":
            flags.append("missing_year")
        if not raw.get("description"):
            flags.append("no_description")
        return flags

    @staticmethod
    def save(records: list[HarmonizedRecord], out_dir: Path) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "harmonized.jsonl"
        with open(out_path, "w") as f:
            for r in records:
                f.write(r.to_jsonl() + "\n")
        logger.info("Saved %d harmonized records → %s", len(records), out_path)
        return out_path
