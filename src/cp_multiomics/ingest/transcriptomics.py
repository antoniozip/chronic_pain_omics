"""Transcriptomics ingester — GEO (Gene Expression Omnibus) via NCBI E-utilities."""

from __future__ import annotations

import logging
import time

from . import register_ingester
from .base import DatasetRecord, OmicsIngester

logger = logging.getLogger(__name__)

GEO_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
GEO_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GEO_FETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GEO_STUDY_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="


@register_ingester("transcriptomics")
class TranscriptomicsIngester(OmicsIngester):
    """Queries GEO for RNA-seq and microarray datasets related to chronic pain."""

    modality = "transcriptomics"

    def _build_query(self, terms: list[str], species_terms: list[str]) -> str:
        pain_clause = " OR ".join(f'"{t}"[Title/Abstract]' for t in terms)
        species_clause = " OR ".join(f'"{s}"[Organism]' for s in species_terms)
        # Restrict to series (GSE) with expression data
        return (
            f"({pain_clause}) AND ({species_clause})"
            ' AND ("Expression profiling by high throughput sequencing"[DataSet Type]'
            ' OR "Expression profiling by array"[DataSet Type])'
            ' AND "gse"[Filter]'
        )

    def _fetch_page(self, query: str, offset: int, page_size: int) -> list[DatasetRecord]:
        # Step 1: search for GEO accession IDs
        search_params = {
            "db": "gds",
            "term": query,
            "retstart": offset,
            "retmax": page_size,
            "retmode": "json",
            "usehistory": "n",
        }
        try:
            search_data = self._get(GEO_SEARCH_URL, params=search_params)
        except Exception as exc:
            logger.error("GEO search failed: %s", exc)
            return []

        ids = search_data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        time.sleep(0.35)  # NCBI rate limit
        return self._fetch_summaries(ids)

    def _fetch_summaries(self, geo_uids: list[str]) -> list[DatasetRecord]:
        """Fetch ESummary for a batch of GEO UIDs and parse into DatasetRecords."""
        params = {
            "db": "gds",
            "id": ",".join(geo_uids),
            "retmode": "json",
        }
        try:
            data = self._get(GEO_SUMMARY_URL, params=params)
        except Exception as exc:
            logger.error("GEO summary fetch failed: %s", exc)
            return []

        records: list[DatasetRecord] = []
        result = data.get("result", {})

        for uid in geo_uids:
            item = result.get(uid)
            if not item or not isinstance(item, dict):
                continue
            accession = item.get("accession", "")
            if not accession.startswith("GSE"):
                continue  # skip GPL/GDS entries, only want series

            n_samples = int(item.get("n_samples", 0) or 0)
            taxon = item.get("taxon", "")
            species = [s.strip() for s in taxon.split(";")] if taxon else []
            platform_info = item.get("gpl", "") or item.get("platform_organism", "")
            year = str(item.get("pdat", ""))[:4]

            records.append(
                DatasetRecord(
                    modality=self.modality,
                    source="GEO",
                    accession=accession,
                    title=item.get("title", ""),
                    description=item.get("summary", ""),
                    species=species,
                    n_samples=n_samples,
                    platform=platform_info,
                    year=year,
                    url=f"{GEO_STUDY_BASE}{accession}",
                    keywords=item.get("subtype", "").split(";") if item.get("subtype") else [],
                )
            )

        return records
