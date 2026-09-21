"""Abstract base class and shared data model for all omics ingesters."""

from __future__ import annotations

import json
import logging
import time
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT = 30
DEFAULT_RETRY_DELAY = 2.0


@dataclass
class DatasetRecord:
    """Normalized metadata record for one omics dataset, regardless of source."""
    modality: str
    source: str        # GEO | GWAS_CATALOG | PRIDE | METABOLIGHTS
    accession: str     # GSE12345, GCST000001, PXD000001, MTBLS123
    title: str
    description: str
    species: list[str]
    n_samples: int
    platform: str      # RNA-seq, microarray, LCMS, etc.
    year: str
    url: str
    keywords: list[str] = field(default_factory=list)
    retrieved_at: str = field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self))


class OmicsIngester(ABC):
    """Base class for all per-modality omics dataset ingesters."""

    modality: str = ""

    def search(self, query: str, filters: dict[str, Any], max_results: int) -> list[DatasetRecord]:
        """Search the repository and return normalized DatasetRecords."""
        raise NotImplementedError

    @abstractmethod
    def _build_query(self, terms: list[str], species_terms: list[str]) -> str:
        """Compose the repository-specific query string."""

    @abstractmethod
    def _fetch_page(self, query: str, offset: int, page_size: int) -> list[DatasetRecord]:
        """Fetch one page of results; must be implemented per source."""

    def search_all(
        self,
        terms: list[str],
        species_terms: list[str],
        filters: dict[str, Any],
        hard_limit: int | None = None,
        page_size: int = 100,
        allow_empty: bool = False,
    ) -> list[DatasetRecord]:
        """Paginate until the source is exhausted.

        Args:
            terms: Topic terms for the repository query.
            species_terms: Species terms for the repository query.
            filters: Post-retrieval filters (applied by the caller).
            hard_limit: Optional safety cap. Reaching it logs a WARNING, because
                a truncated result set is not the same as an exhausted one.
            page_size: Records per request.
            allow_empty: Permit a zero-record result instead of raising.

        Returns:
            All DatasetRecords the source served.

        Raises:
            RuntimeError: If the source returned no records and allow_empty is False.
                A silent zero is how the Europe PMC and PRIDE keyword defects went
                unnoticed for months.
        """
        query = self._build_query(terms, species_terms)
        logger.info("[%s] Query: %s", self.modality, query[:200])

        records: list[DatasetRecord] = []
        offset = 0

        while True:
            if hard_limit is not None and len(records) >= hard_limit:
                logger.warning(
                    "[%s] hard_limit reached (%d): results are TRUNCATED, not exhausted",
                    self.modality, hard_limit,
                )
                break

            want = page_size
            if hard_limit is not None:
                want = min(page_size, hard_limit - len(records))

            batch = self._fetch_page(query, offset, want)
            if not batch:
                break

            records.extend(batch)
            logger.info("[%s] Fetched %d records", self.modality, len(records))
            offset += len(batch)

            if len(batch) < want:
                break
            time.sleep(0.3)

        if not records and not allow_empty:
            raise RuntimeError(
                f"[{self.modality}] search returned 0 records for query: {query[:200]!r}. "
                "This is almost always a broken endpoint or an unsupported query syntax, "
                "not a genuinely empty repository. Pass allow_empty=True to override."
            )

        return records

    # ------------------------------------------------------------------
    # HTTP helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _get(url: str, params: dict | None = None, timeout: int = DEFAULT_TIMEOUT) -> Any:
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        return resp.json()

    @staticmethod
    def save(records: list[DatasetRecord], out_dir: Path) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "datasets.jsonl"
        with open(out_path, "w") as f:
            for r in records:
                f.write(r.to_jsonl() + "\n")
        logger.info("Saved %d records → %s", len(records), out_path)
        return out_path
