"""Genomics ingester — GWAS Catalog REST API."""

from __future__ import annotations

import logging
from typing import Any

from . import register_ingester
from .base import DatasetRecord, OmicsIngester

logger = logging.getLogger(__name__)

GWAS_API_BASE = "https://www.ebi.ac.uk/gwas/rest/api"
GWAS_STUDIES_URL = f"{GWAS_API_BASE}/studies/search/findByDiseaseTrait"
GWAS_FREE_TEXT_URL = f"{GWAS_API_BASE}/studies/search/findByDiseaseTrait"
GWAS_STUDY_UI = "https://www.ebi.ac.uk/gwas/studies/"

# GWAS Catalog EFO terms for pain-related traits
PAIN_EFO_TERMS = [
    "chronic pain",
    "neuropathic pain",
    "back pain",
    "fibromyalgia",
    "migraine",
    "pain measurement",
    "musculoskeletal pain",
]


@register_ingester("genomics")
class GenomicsIngester(OmicsIngester):
    """Queries GWAS Catalog for GWAS studies on chronic pain traits."""

    modality = "genomics"

    def _build_query(self, terms: list[str], species_terms: list[str]) -> str:
        # GWAS Catalog uses free-text trait search; combine user terms with known EFOs
        all_terms = list(set(terms + PAIN_EFO_TERMS))
        return "|".join(all_terms)  # we search each term separately

    def _fetch_page(self, query: str, offset: int, page_size: int) -> list[DatasetRecord]:
        # GWAS Catalog paginates differently (page number, not offset)
        page = offset // page_size
        terms = [t.strip() for t in query.split("|") if t.strip()]
        # Search the first term only when paginating; multi-term handled in search_all override
        term = terms[0] if terms else "chronic pain"

        params = {
            "diseaseTrait": term,
            "page": page,
            "size": page_size,
        }
        try:
            data = self._get(GWAS_STUDIES_URL, params=params)
        except Exception as exc:
            logger.error("GWAS Catalog fetch failed for term '%s': %s", term, exc)
            return []

        return self._parse_response(data)

    def search_all(
        self,
        terms: list[str],
        species_terms: list[str],
        filters: dict[str, Any],
        hard_limit: int | None = None,
        page_size: int = 100,
        allow_empty: bool = False,
    ) -> list[DatasetRecord]:
        """Override to iterate over each pain-related EFO term."""
        all_records: list[DatasetRecord] = []
        seen_accessions: set[str] = set()

        search_terms = list(set(terms + PAIN_EFO_TERMS))

        for term in search_terms:
            if hard_limit is not None and len(all_records) >= hard_limit:
                break
            page = 0
            while hard_limit is None or len(all_records) < hard_limit:
                params = {"diseaseTrait": term, "page": page, "size": page_size}
                try:
                    data = self._get(GWAS_STUDIES_URL, params=params)
                except Exception as exc:
                    logger.warning("GWAS term '%s' page %d failed: %s", term, page, exc)
                    break

                batch = self._parse_response(data)
                if not batch:
                    break

                for r in batch:
                    if r.accession not in seen_accessions:
                        seen_accessions.add(r.accession)
                        all_records.append(r)

                embedded = data.get("_embedded", {})
                studies = embedded.get("studies", [])
                if len(studies) < page_size:
                    break
                page += 1

            logger.info("[genomics] Term '%s': %d unique studies so far", term, len(all_records))

        if hard_limit is not None and len(all_records) >= hard_limit:
            logger.warning(
                "[genomics] hard_limit reached (%d): results are TRUNCATED, not exhausted",
                hard_limit,
            )

        if not all_records and not allow_empty:
            raise RuntimeError(
                f"[{self.modality}] search returned 0 records for terms {search_terms!r}. "
                "This is almost always a broken endpoint or an unsupported query syntax, "
                "not a genuinely empty repository. Pass allow_empty=True to override."
            )

        return all_records if hard_limit is None else all_records[:hard_limit]

    def _parse_response(self, data: dict) -> list[DatasetRecord]:
        embedded = data.get("_embedded", {})
        studies = embedded.get("studies", [])
        records: list[DatasetRecord] = []

        for s in studies:
            accession = s.get("accessionId", "")
            if not accession:
                continue
            trait = s.get("diseaseTrait", {}).get("trait", "")
            pub = s.get("publicationInfo", {})
            year = str(pub.get("publicationDate", ""))[:4]
            # Best-effort sample count from the descriptor string
            sample_count = self._parse_sample_count(s.get("initialSampleSize", ""))

            records.append(
                DatasetRecord(
                    modality=self.modality,
                    source="GWAS_CATALOG",
                    accession=accession,
                    title=pub.get("title", ""),
                    description=trait,
                    species=["Homo sapiens"],  # GWAS Catalog is human-only
                    n_samples=sample_count,
                    platform="GWAS",
                    year=year,
                    url=f"{GWAS_STUDY_UI}{accession}",
                    keywords=[trait],
                )
            )

        return records

    @staticmethod
    def _parse_sample_count(descriptor: str) -> int:
        """Extract the first integer from a sample descriptor like '4,326 European'."""
        import re
        m = re.search(r"[\d,]+", descriptor or "")
        if not m:
            return 0
        try:
            return int(m.group().replace(",", ""))
        except ValueError:
            return 0
