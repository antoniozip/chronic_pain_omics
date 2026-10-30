"""Ortholog mapping via MyGene.info REST API with local disk cache.

Maps gene symbols from mouse (Mus musculus, taxid 10090) or rat
(Rattus norvegicus, taxid 10116) to human (Homo sapiens, taxid 9606)
HGNC symbols. Uses MyGene.info POST /v3/query + /v3/gene endpoints.

Two-step mapping:
  1. POST /v3/query: animal symbol → homologene data (contains human Entrez ID)
  2. POST /v3/gene: human Entrez IDs → human HGNC symbols

Cache: data/interim/ortholog_cache.json (keyed by "{symbol}:{from_taxid}")
  Versioned (`__meta__`); an unversioned cache is rejected because older
  versions wrote `None` for a transient API failure, which is permanent.
  A `None` is only cached when the service explicitly reports no homologene
  data; request failures are retried and then propagated.
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
#: MyGene rejects `size` at 400 with 400 "Too Many Hits Requested".
#: Keep the batch comfortably below the cap; each batch is one query term set
#: and `size` is set to the batch length so no hit is silently dropped.
BATCH_SIZE = 200
RATE_DELAY = 0.5   # seconds between batches

#: Cache format. Older files are flat `{key: record|None}` maps, and a None in
#: one of them may be a transient API failure that was written as "no ortholog"
#: forever. They are rejected and rebuilt instead of trusted.
CACHE_VERSION = 2
CACHE_META_KEY = "__meta__"

#: A failed request is retried this many times before propagating. The point is
#: that it never becomes a cached None: a timeout must not remove a gene from
#: every future cross-species run.
RETRY_ATTEMPTS = 3

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

        Returns dict keyed by input symbol; value is OrthologRecord, or None
        when no ortholog was found or the service did not answer (the latter is
        not cached and is retried on the next run).
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

    def _post_json(self, url: str, payload: dict, what: str):
        """POST and parse JSON, retrying a bounded number of times.

        A request or parse failure propagates as RuntimeError. Returning an
        empty result would be cached as "no ortholog", converting a transient
        outage into a permanent hole in the cross-species comparison.
        """
        last: Exception | None = None
        for attempt in range(RETRY_ATTEMPTS):
            try:
                resp = requests.post(url, data=payload, timeout=60)
                resp.raise_for_status()
                return resp.json()
            except Exception as exc:
                last = exc
                if attempt + 1 < RETRY_ATTEMPTS:
                    delay = RATE_DELAY * (2 ** attempt)
                    logger.warning(
                        "%s failed (attempt %d/%d): %s; retrying in %.1fs",
                        what, attempt + 1, RETRY_ATTEMPTS, exc, delay,
                    )
                    time.sleep(delay)
        raise RuntimeError(
            f"{what} failed after {RETRY_ATTEMPTS} attempts: {last}"
        ) from last

    def _query_mygene(
        self, symbols: list[str], from_taxid: int
    ) -> dict[str, OrthologRecord | None]:
        # Step 1: symbol → homologene (contains human Entrez IDs)
        payload = {
            "q": ",".join(symbols),
            "scopes": "symbol",
            "species": str(from_taxid),
            "fields": "symbol,homologene",
            "size": str(len(symbols)),
        }
        hits = self._post_json(MYGENE_QUERY_URL, payload, "MyGene.info query")
        if not isinstance(hits, list):
            raise RuntimeError(
                "MyGene.info query returned "
                f"{type(hits).__name__}, expected a list of hits"
            )

        # Extract human Entrez IDs from homologene data. Only symbols the
        # service actually answered about may be cached as absent; a symbol
        # missing from the response is left uncached so it is retried.
        addressed: set[str] = set()
        sym_to_entrez: dict[str, int] = {}
        for hit in hits:
            if not isinstance(hit, dict):
                continue
            query_sym = hit.get("query", "")
            if not query_sym:
                continue
            addressed.add(query_sym)
            if hit.get("notfound"):
                continue
            genes = (hit.get("homologene") or {}).get("genes", []) or []
            for taxid_val, entrez_id in genes:
                if taxid_val == TAXID["Homo sapiens"]:
                    sym_to_entrez[query_sym] = entrez_id
                    break

        result: dict[str, OrthologRecord | None] = {s: None for s in symbols
                                                   if s in addressed}
        unanswered = [s for s in symbols if s not in addressed]
        if unanswered:
            logger.warning(
                "%d of %d symbols were not answered by MyGene.info; leaving "
                "them uncached so they are retried (for example: %s)",
                len(unanswered), len(symbols), ", ".join(unanswered[:5]),
            )
        if not sym_to_entrez:
            return result

        # Step 2: resolve human Entrez IDs → HGNC symbols
        entrez_ids = sorted(set(sym_to_entrez.values()))
        gene_hits = self._post_json(
            MYGENE_GENE_URL,
            {"ids": ",".join(str(e) for e in entrez_ids), "fields": "symbol,taxid"},
            "MyGene.info gene lookup",
        )
        if not isinstance(gene_hits, list):
            raise RuntimeError(
                "MyGene.info gene lookup returned "
                f"{type(gene_hits).__name__}, expected a list of hits"
            )

        entrez_to_symbol: dict[int, str] = {}
        for g in gene_hits:
            if (isinstance(g, dict)
                    and g.get("taxid") == TAXID["Homo sapiens"]
                    and "symbol" in g):
                entrez_to_symbol[int(g["_id"])] = g["symbol"]

        species_name = next((s for s, t in TAXID.items() if t == from_taxid), str(from_taxid))
        for sym, entrez in sym_to_entrez.items():
            human_sym = entrez_to_symbol.get(entrez, "")
            if not human_sym:
                # The service answered and carried no human symbol: that is a
                # confirmed absence, not a failure.
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
            meta = raw.get(CACHE_META_KEY) if isinstance(raw, dict) else None
            version = meta.get("version") if isinstance(meta, dict) else None
            if version != CACHE_VERSION:
                logger.warning(
                    "Ignoring ortholog cache %s: format version %r (expected %d). "
                    "Rebuilding; a legacy cache can hold None entries written by "
                    "a failed fetch, which would be permanent.",
                    self._cache_path, version, CACHE_VERSION,
                )
                return
            for key, val in raw.items():
                if key == CACHE_META_KEY:
                    continue
                try:
                    self._cache[key] = None if val is None else OrthologRecord(**val)
                except (TypeError, ValueError) as exc:
                    # One malformed record must not truncate the load: the next
                    # save would then replace the file with that prefix.
                    logger.warning(
                        "Dropping malformed ortholog cache entry %s: %s", key, exc)
            logger.debug("Ortholog cache loaded: %d entries", len(self._cache))
        except Exception as exc:
            logger.warning("Could not load ortholog cache: %s", exc)

    def _save_cache(self) -> None:
        self._cache_path.parent.mkdir(parents=True, exist_ok=True)
        serialized: dict[str, object] = {
            CACHE_META_KEY: {"version": CACHE_VERSION},
        }
        serialized.update({
            k: (None if v is None else v.__dict__)
            for k, v in self._cache.items()
        })
        tmp = self._cache_path.with_suffix(self._cache_path.suffix + ".tmp")
        with open(tmp, "w") as f:
            json.dump(serialized, f, indent=2)
        tmp.replace(self._cache_path)
