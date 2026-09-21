"""Ortholog mapping via MyGene.info REST API with local disk cache.

Maps gene symbols from mouse (Mus musculus, taxid 10090) or rat
(Rattus norvegicus, taxid 10116) to human (Homo sapiens, taxid 9606)
HGNC symbols. Uses MyGene.info POST /v3/query + /v3/gene endpoints.

Two-step mapping:
  1. POST /v3/query: animal symbol → homologene data (contains human Entrez ID)
  2. POST /v3/gene: human Entrez IDs → human HGNC symbols

Cache: data/interim/ortholog_cache.json (keyed by "{symbol}:{from_taxid}")
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

MYGENE_QUERY_URL = "https://mygene.info/v3/query"
MYGENE_GENE_URL  = "https://mygene.info/v3/gene"
BATCH_SIZE = 500
RATE_DELAY = 0.5   # seconds between batches

TAXID = {
    "Homo sapiens":       9606,
    "Mus musculus":       10090,
    "Rattus norvegicus":  10116,
}

SUPPORTED_ANIMAL_SPECIES = {"Mus musculus", "Rattus norvegicus"}


@dataclass
class OrthologRecord:
    source_symbol: str
    source_species: str
    human_symbol: str
    human_entrez: int | None
    confidence: str        # "high" | "medium" | "low" (based on n sources)
    sources: list[str] = field(default_factory=list)


class OrthologMapper:
    """Maps animal gene symbols to human HGNC symbols via MyGene.info.

    Args:
        cache_path: Path to JSON cache file (created on first use).
    """

    def __init__(self, cache_path: Path) -> None:
        self._cache_path = cache_path
        self._cache: dict[str, OrthologRecord | None] = {}
        self._load_cache()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def map_to_human(
        self,
        symbols: list[str],
        from_species: str,
    ) -> dict[str, OrthologRecord | None]:
        """Map a list of animal gene symbols to human orthologs.

        Returns dict keyed by input symbol; value is OrthologRecord or None
        if no ortholog was found.
        """
        if from_species not in SUPPORTED_ANIMAL_SPECIES:
            raise ValueError(
                f"Unsupported species '{from_species}'. "
                f"Supported: {sorted(SUPPORTED_ANIMAL_SPECIES)}"
            )

        taxid = TAXID[from_species]
        cache_keys = {s: f"{s}:{taxid}" for s in symbols}
        missing = [s for s, k in cache_keys.items() if k not in self._cache]

        if missing:
            logger.info(
                "Querying MyGene.info for %d unmapped symbols (%s → human)",
                len(missing), from_species,
            )
            self._fetch_and_cache(missing, taxid)
            self._save_cache()

        return {s: self._cache.get(cache_keys[s]) for s in symbols}

    # ------------------------------------------------------------------
    # MyGene.info fetch
    # ------------------------------------------------------------------

    def _fetch_and_cache(self, symbols: list[str], from_taxid: int) -> None:
        for i in range(0, len(symbols), BATCH_SIZE):
            batch = symbols[i : i + BATCH_SIZE]
            results = self._query_mygene(batch, from_taxid)
            for symbol, record in results.items():
                cache_key = f"{symbol}:{from_taxid}"
                self._cache[cache_key] = record
            if i + BATCH_SIZE < len(symbols):
                time.sleep(RATE_DELAY)

    def _query_mygene(
        self, symbols: list[str], from_taxid: int
    ) -> dict[str, OrthologRecord | None]:
        # Step 1: symbol → homologene (contains human Entrez IDs)
        payload = {
            "q": ",".join(symbols),
            "scopes": "symbol",
            "species": str(from_taxid),
            "fields": "symbol,homologene",
            "size": str(BATCH_SIZE),
        }
        try:
            resp = requests.post(MYGENE_QUERY_URL, data=payload, timeout=60)
            resp.raise_for_status()
            hits = resp.json()
        except Exception as exc:
            logger.error("MyGene.info query failed: %s", exc)
            return {s: None for s in symbols}

        # Extract human Entrez IDs from homologene data
        sym_to_entrez: dict[str, int] = {}
        for hit in hits:
            query_sym = hit.get("query", "")
            if hit.get("notfound") or "homologene" not in hit:
                continue
            genes = hit["homologene"].get("genes", [])
            for taxid_val, entrez_id in genes:
                if taxid_val == TAXID["Homo sapiens"]:
                    sym_to_entrez[query_sym] = entrez_id
                    break

        if not sym_to_entrez:
            return {s: None for s in symbols}

        # Step 2: resolve human Entrez IDs → HGNC symbols
        entrez_ids = list(set(sym_to_entrez.values()))
        entrez_to_symbol: dict[int, str] = {}
        try:
            gene_resp = requests.post(
                MYGENE_GENE_URL,
                data={"ids": ",".join(str(e) for e in entrez_ids), "fields": "symbol,taxid"},
                timeout=60,
            )
            gene_resp.raise_for_status()
            for g in gene_resp.json():
                if g.get("taxid") == TAXID["Homo sapiens"] and "symbol" in g:
                    entrez_to_symbol[int(g["_id"])] = g["symbol"]
        except Exception as exc:
            logger.error("MyGene.info gene lookup failed: %s", exc)

        species_name = next((s for s, t in TAXID.items() if t == from_taxid), str(from_taxid))
        result: dict[str, OrthologRecord | None] = {s: None for s in symbols}

        for sym, entrez in sym_to_entrez.items():
            human_sym = entrez_to_symbol.get(entrez, "")
            if not human_sym:
                continue
            result[sym] = OrthologRecord(
                source_symbol=sym,
                source_species=species_name,
                human_symbol=human_sym,
                human_entrez=entrez,
                confidence="high",
                sources=["homologene"],
            )

        return result

    # ------------------------------------------------------------------
    # Cache persistence
    # ------------------------------------------------------------------

    def _load_cache(self) -> None:
        if not self._cache_path.exists():
            return
        try:
            with open(self._cache_path) as f:
                raw = json.load(f)
            for key, val in raw.items():
                if val is None:
                    self._cache[key] = None
                else:
                    self._cache[key] = OrthologRecord(**val)
            logger.debug("Ortholog cache loaded: %d entries", len(self._cache))
        except Exception as exc:
            logger.warning("Could not load ortholog cache: %s", exc)

    def _save_cache(self) -> None:
        self._cache_path.parent.mkdir(parents=True, exist_ok=True)
        serialized = {
            k: (None if v is None else v.__dict__)
            for k, v in self._cache.items()
        }
        with open(self._cache_path, "w") as f:
            json.dump(serialized, f, indent=2)
