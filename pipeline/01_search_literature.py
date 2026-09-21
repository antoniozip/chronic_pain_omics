"""Literature search pipeline step.

Queries PubMed and Europe PMC using the config in conf/search/chronic_pain.yaml.
Results are deduplicated by PMID/DOI and written to data/raw/literature/.
All search decisions are logged to literature/prisma/search_log.csv (PRISMA tracking).

Usage:
    uv run python pipeline/01_search_literature.py
    uv run python pipeline/01_search_literature.py --query-config conf/search/my_query.yaml
    uv run python pipeline/01_search_literature.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import requests
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
RAW_LITERATURE_DIR = REPO_ROOT / "data" / "raw" / "literature"
PRISMA_LOG = REPO_ROOT / "literature" / "prisma" / "search_log.csv"

PUBMED_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_FETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EUROPE_PMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

@dataclass
class SearchRecord:
    """One database hit, normalized across sources."""
    source: str          # "pubmed" | "europe_pmc"
    uid: str             # PMID or Europe PMC ID
    doi: str
    title: str
    authors: str         # semicolon-joined
    journal: str
    year: str
    abstract: str
    query_label: str     # which sub-query produced this hit
    retrieved_at: str    # ISO-8601

    def dedup_key(self) -> str:
        """Stable key for deduplication: prefer DOI, fall back to uid+source."""
        if self.doi:
            return f"doi:{self.doi.lower().strip()}"
        return f"{self.source}:{self.uid}"


@dataclass
class PrismaLog:
    """Running PRISMA flow totals for one search run."""
    run_id: str
    config_file: str
    databases: list[str]
    total_hits: int = 0
    duplicates_removed: int = 0
    after_dedup: int = 0
    timestamp: str = field(default_factory=lambda: datetime.now(UTC).isoformat())


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(config_path: Path) -> dict[str, Any]:
    with open(config_path) as f:
        return yaml.safe_load(f)


def build_topic_query(cfg: dict[str, Any]) -> str:
    """Build the source-agnostic topic/omics/species clause (no date filter)."""
    terms = cfg["query"]["terms"]
    omics_parts = [v for v in cfg["query"]["omics"].values()]
    species_parts = list(cfg["query"]["species"].values())

    topic = " AND ".join(f"({t})" for t in terms)
    omics = " OR ".join(f"({o})" for o in omics_parts)
    species = " OR ".join(f"({s})" for s in species_parts)

    return f"({topic}) AND ({omics}) AND ({species})"


def build_query(cfg: dict[str, Any], source: str) -> str:
    """Build a boolean query for one bibliographic source.

    The date filter is source-specific. PubMed uses the [PDAT] field tag; Europe
    PMC cannot parse it and silently returns hitCount=0 for any query containing
    it, which is why Europe PMC contributed zero records to every run before
    2026-07-10.

    Args:
        cfg: Parsed search config.
        source: "pubmed" or "europe_pmc".

    Returns:
        The full query string for that source.

    Raises:
        ValueError: If source is not a supported database.
    """
    if source not in ("pubmed", "europe_pmc"):
        raise ValueError(f"unknown source: {source!r}")

    query = build_topic_query(cfg)

    date_start = cfg.get("filters", {}).get("date_range", {}).get("start", "")
    if not date_start:
        return query

    date_end = cfg.get("filters", {}).get("date_range", {}).get("end", "") or "2100-12-31"

    if source == "pubmed":
        return f'{query} AND ("{date_start[:4]}"[PDAT]:"{str(date_end)[:4]}"[PDAT])'

    return f"{query} AND (FIRST_PDATE:[{date_start} TO {date_end}])"


def build_pubmed_query(cfg: dict[str, Any]) -> str:
    """Backwards-compatible alias for build_query(cfg, "pubmed")."""
    return build_query(cfg, "pubmed")


# ---------------------------------------------------------------------------
# PubMed retrieval
# ---------------------------------------------------------------------------

def search_pubmed(query: str, max_results: int, dry_run: bool) -> list[str]:
    """Return list of PMIDs matching query."""
    if dry_run:
        logger.info("[dry-run] Would search PubMed: %s", query[:120])
        return []

    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "usehistory": "y",
        "retmode": "json",
    }
    resp = requests.get(PUBMED_SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    ids = data["esearchresult"]["idlist"]
    logger.info("PubMed: %d hits", len(ids))
    return ids


def fetch_pubmed_records(
    pmids: list[str], query_label: str, batch_size: int = 200
) -> list[SearchRecord]:
    """Fetch PubMed abstracts in batches and parse into SearchRecords."""
    records: list[SearchRecord] = []
    now = datetime.now(UTC).isoformat()

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "abstract",
            "retmode": "xml",
        }
        resp = requests.get(PUBMED_FETCH_URL, params=params, timeout=60)
        resp.raise_for_status()
        records.extend(_parse_pubmed_xml(resp.text, query_label, now))
        time.sleep(0.4)  # respect NCBI rate limit (3 req/s without API key)

    return records


def _parse_pubmed_xml(xml_text: str, query_label: str, retrieved_at: str) -> list[SearchRecord]:
    """Minimal XML parse without lxml dependency."""
    from xml.etree import ElementTree as ET

    root = ET.fromstring(xml_text)
    records: list[SearchRecord] = []

    for article in root.iter("PubmedArticle"):
        try:
            pmid = article.findtext(".//PMID", "")
            title = article.findtext(".//ArticleTitle", "")
            abstract = " ".join(
                t.text or "" for t in article.findall(".//AbstractText")
            )
            journal = article.findtext(".//Journal/Title", "")
            year = article.findtext(".//PubDate/Year", "") or article.findtext(
                ".//PubDate/MedlineDate", ""
            )[:4]
            authors = "; ".join(
                f"{a.findtext('LastName', '')} {a.findtext('Initials', '')}".strip()
                for a in article.findall(".//Author")
            )
            doi = next(
                (
                    id_el.text or ""
                    for id_el in article.findall(".//ArticleId")
                    if id_el.get("IdType") == "doi"
                ),
                "",
            )
            records.append(
                SearchRecord(
                    source="pubmed",
                    uid=pmid,
                    doi=doi,
                    title=title,
                    authors=authors,
                    journal=journal,
                    year=year,
                    abstract=abstract,
                    query_label=query_label,
                    retrieved_at=retrieved_at,
                )
            )
        except Exception as exc:
            logger.warning("Failed to parse PubMed article: %s", exc)

    return records


# ---------------------------------------------------------------------------
# Europe PMC retrieval
# ---------------------------------------------------------------------------

def search_europe_pmc(
    query: str, max_results: int, query_label: str, dry_run: bool
) -> list[SearchRecord]:
    """Page through Europe PMC REST API and return SearchRecords."""
    if dry_run:
        logger.info("[dry-run] Would search Europe PMC: %s", query[:120])
        return []

    records: list[SearchRecord] = []
    now = datetime.now(UTC).isoformat()
    cursor = "*"
    page_size = 1000

    while len(records) < max_results:
        params = {
            "query": query,
            "resultType": "core",
            "format": "json",
            "pageSize": min(page_size, max_results - len(records)),
            "cursorMark": cursor,
        }
        resp = requests.get(EUROPE_PMC_URL, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        hits = data.get("resultList", {}).get("result", [])
        if not hits:
            break

        for h in hits:
            records.append(
                SearchRecord(
                    source="europe_pmc",
                    uid=h.get("id", ""),
                    doi=h.get("doi", ""),
                    title=h.get("title", ""),
                    authors=h.get("authorString", ""),
                    journal=h.get("journalTitle", ""),
                    year=str(h.get("pubYear", "")),
                    abstract=h.get("abstractText", ""),
                    query_label=query_label,
                    retrieved_at=now,
                )
            )

        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        time.sleep(0.2)

    logger.info("Europe PMC: %d records", len(records))
    return records


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate(records: list[SearchRecord]) -> tuple[list[SearchRecord], int]:
    """Remove duplicates by DOI (preferred) or source+uid. Returns (unique, n_removed)."""
    seen: set[str] = set()
    unique: list[SearchRecord] = []
    for r in records:
        key = r.dedup_key()
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique, len(records) - len(unique)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_records(records: list[SearchRecord], out_dir: Path, run_id: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{run_id}_hits.jsonl"
    with open(out_path, "w") as f:
        for r in records:
            f.write(json.dumps(asdict(r)) + "\n")
    logger.info("Saved %d records → %s", len(records), out_path)
    return out_path


def append_prisma_log(prisma: PrismaLog) -> None:
    PRISMA_LOG.parent.mkdir(parents=True, exist_ok=True)
    write_header = not PRISMA_LOG.exists()
    with open(PRISMA_LOG, "a", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "run_id", "config_file", "databases",
                "total_hits", "duplicates_removed", "after_dedup", "timestamp",
            ],
        )
        if write_header:
            writer.writeheader()
        row = asdict(prisma)
        row["databases"] = ";".join(row["databases"])
        writer.writerow(row)
    logger.info("PRISMA log updated → %s", PRISMA_LOG)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--query-config",
        type=Path,
        default=REPO_ROOT / "conf" / "search" / "chronic_pain.yaml",
        help="Path to search YAML config",
    )
    p.add_argument("--dry-run", action="store_true", help="Print queries without hitting APIs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.query_config)

    run_id = hashlib.sha1(
        (str(args.query_config) + datetime.now(UTC).isoformat()).encode()
    ).hexdigest()[:10]

    databases: list[str] = cfg.get("databases", ["pubmed"])
    max_results: int = cfg.get("output", {}).get("max_results_per_query", 5000)
    query_label = args.query_config.stem

    all_records: list[SearchRecord] = []

    if "pubmed" in databases:
        query = build_query(cfg, "pubmed")
        logger.info("PubMed query: %s", query[:200])
        pmids = search_pubmed(query, max_results, dry_run=args.dry_run)
        if pmids:
            all_records.extend(fetch_pubmed_records(pmids, query_label=query_label))

    if "europe_pmc" in databases:
        query = build_query(cfg, "europe_pmc")
        logger.info("Europe PMC query: %s", query[:200])
        all_records.extend(
            search_europe_pmc(query, max_results, query_label=query_label, dry_run=args.dry_run)
        )

    total_hits = len(all_records)
    unique_records, n_dupes = deduplicate(all_records)

    prisma = PrismaLog(
        run_id=run_id,
        config_file=str(args.query_config),
        databases=databases,
        total_hits=total_hits,
        duplicates_removed=n_dupes,
        after_dedup=len(unique_records),
    )

    logger.info(
        "Total: %d | Duplicates removed: %d | After dedup: %d",
        total_hits, n_dupes, len(unique_records),
    )

    if not args.dry_run:
        out_dir = Path(cfg.get("output", {}).get("raw_hits_dir", RAW_LITERATURE_DIR))
        save_records(unique_records, out_dir, run_id)
        append_prisma_log(prisma)
    else:
        logger.info("[dry-run] Skipped writing output files.")


if __name__ == "__main__":
    main()
