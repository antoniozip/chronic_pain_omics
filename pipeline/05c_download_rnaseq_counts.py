"""Download RNA-seq count matrices from GEO supplementary files.

For studies that have a _groups.csv annotation but no per-study effects file,
this script:
  1. Queries the GEO accession text API to list supplementary files.
  2. Identifies likely count matrix files (patterns: "count", "matrix",
     "raw", "expression", "gene_exp").
  3. Downloads them to data/raw/geo_cache/{accession}/.

Skips studies that already have supplementary files downloaded.

Usage:
  python pipeline/05c_download_rnaseq_counts.py
  python pipeline/05c_download_rnaseq_counts.py --accessions GSE102721 GSE111216
  python pipeline/05c_download_rnaseq_counts.py --dry-run
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import time
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = REPO_ROOT / "data" / "raw" / "geo_cache"

GEO_TEXT_API = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# File patterns that suggest a count / expression matrix (case-insensitive)
COUNT_PATTERNS = re.compile(
    r"(count|raw.count|gene.count|expression.matrix|gene.exp|read.count"
    r"|feature.count|htseq|deseq|edger|rna.seq|counts.matrix"
    r"|count_matrix|raw_counts|gene_counts|read_counts"
    r"|tpm|rpkm|fpkm|cpm|normalized.count"
    r"|expression.gene|gene.expression|read.per.feature|reads.per"
    r"|matrix\.mtx|sparse.matrix"
    r"|allsamples.*count|combined.*count|bulk.*count"
    r"|\.count\.|_counts\."
    r"|expression_data|exprs|expr_mat)",
    re.IGNORECASE,
)

# Extensions we want (add .mtx for 10x sparse matrices)
WANTED_EXT = {".gz", ".zip", ".txt", ".csv", ".tsv", ".xlsx", ".xls", ".tar", ".mtx", ".rds"}

# Files to skip (large alignment files, index files, etc.)
SKIP_PATTERNS = re.compile(
    r"\.(bam|bai|cram|crai|sam|sra|fastq|fq|bigwig|bw|bed|gtf|gff|vcf)(\.gz)?$",
    re.IGNORECASE,
)

# Max file size to download (200 MB) — avoid huge tarball downloads
MAX_BYTES = 200 * 1024 * 1024

RATE_DELAY = 1.0  # seconds between GEO API calls


# ---------------------------------------------------------------------------
# GEO API helpers
# ---------------------------------------------------------------------------

def geo_text(accession: str) -> str | None:
    """Fetch the GEO accession text record."""
    url = f"{GEO_TEXT_API}?acc={accession}&targ=self&form=text&view=brief"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.text
    except Exception as exc:
        logger.warning("GEO text API failed for %s: %s", accession, exc)
        return None


def parse_suppl_files(text: str) -> list[str]:
    """Extract supplementary file FTP URLs from GEO text record."""
    urls: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if "!Series_supplementary_file" in line or "!Sample_supplementary_file" in line:
            # Format: !Series_supplementary_file = ftp://...
            parts = line.split("=", 1)
            if len(parts) == 2:
                url = parts[1].strip()
                if url.startswith("ftp://") or url.startswith("https://"):
                    urls.append(url)
    return urls


def ftp_to_https(url: str) -> str:
    """Convert ftp:// GEO URL to https:// equivalent."""
    return url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")


def is_count_file(url: str) -> bool:
    """Return True if the URL looks like a count matrix file."""
    name = Path(url).name
    if SKIP_PATTERNS.search(name):
        return False
    # Suffix check (accounting for .txt.gz etc.)
    suffixes = "".join(Path(name).suffixes).lower()
    if not any(ext in suffixes for ext in WANTED_EXT):
        return False
    return bool(COUNT_PATTERNS.search(name))


def file_size(url: str) -> int | None:
    """Return remote file size from HEAD request, or None."""
    try:
        resp = requests.head(ftp_to_https(url), timeout=10, allow_redirects=True)
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl else None
    except Exception:
        return None


def download_file(url: str, dest: Path) -> bool:
    """Download a file with progress logging. Returns True on success."""
    https_url = ftp_to_https(url)
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(https_url, stream=True, timeout=120) as resp:
            resp.raise_for_status()
            downloaded = 0
            with open(dest, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=65536):
                    fh.write(chunk)
                    downloaded += len(chunk)
            logger.info("    Downloaded %s (%.1f MB)", dest.name, downloaded / 1e6)
        return True
    except Exception as exc:
        logger.warning("    Download failed for %s: %s", dest.name, exc)
        if dest.exists():
            dest.unlink()
        return False


# ---------------------------------------------------------------------------
# Study selection helpers
# ---------------------------------------------------------------------------

def _get_platform(accession: str, cache_dir: Path) -> str | None:
    for path in sorted(cache_dir.glob(f"{accession}*_series_matrix.txt.gz")):
        try:
            with gzip.open(path, "rt", errors="replace") as fh:
                for line in fh:
                    if line.startswith("!series_matrix_table_begin"):
                        break
                    if line.startswith("!Sample_platform_id"):
                        return line.split("\t")[1].strip().strip('"')
        except OSError:
            pass
    return None


RNA_SEQ_GPLS = {
    "GPL10287", "GPL10669", "GPL11154", "GPL13112", "GPL14844", "GPL16417",
    "GPL17021", "GPL18573", "GPL18694", "GPL19057", "GPL20084", "GPL20301",
    "GPL21103", "GPL21163", "GPL22396", "GPL23479", "GPL23945", "GPL24247",
    "GPL24676", "GPL24688", "GPL24782", "GPL25915", "GPL25947", "GPL27943",
    "GPL28457", "GPL32027", "GPL33896", "GPL34290",
}


def find_rna_seq_candidates(cache_dir: Path, per_study_dir: Path) -> list[str]:
    """Return accessions that are annotated, RNA-seq, and lack per-study effects."""
    annotated = {p.name.replace("_groups.csv", "") for p in cache_dir.glob("*_groups.csv")}
    processed = {
        p.name.replace("_effects.csv", "").replace("_effects_normalized", "")
        for p in per_study_dir.glob("*_effects*.csv")
        if "normalized_all" not in p.name
    }
    skipped = annotated - processed
    candidates = []
    for acc in sorted(skipped):
        gpl = _get_platform(acc, cache_dir)
        if gpl in RNA_SEQ_GPLS:
            candidates.append(acc)
    return candidates


def already_downloaded(accession: str, cache_dir: Path) -> bool:
    """Return True if a count-like file exists in the accession subdirectory."""
    suppl_dir = cache_dir / accession
    if not suppl_dir.exists():
        return False
    for f in suppl_dir.iterdir():
        if COUNT_PATTERNS.search(f.name) and not SKIP_PATTERNS.search(f.name):
            return True
    return False


# ---------------------------------------------------------------------------
# Main download logic per study
# ---------------------------------------------------------------------------

def process_study(accession: str, cache_dir: Path, dry_run: bool = False) -> dict:
    record: dict = {"accession": accession, "status": "unknown",
                    "downloaded": [], "skipped": []}

    if already_downloaded(accession, cache_dir):
        record["status"] = "already_present"
        logger.info("%s: supplementary files already present — skipping", accession)
        return record

    logger.info("%s: querying GEO ...", accession)
    text = geo_text(accession)
    time.sleep(RATE_DELAY)
    if not text:
        record["status"] = "api_error"
        return record

    urls = parse_suppl_files(text)
    if not urls:
        logger.info("  %s: no supplementary files listed", accession)
        record["status"] = "no_suppl"
        return record

    count_urls = [u for u in urls if is_count_file(u)]
    other_urls  = [u for u in urls if not is_count_file(u)]
    logger.info("  %s: %d suppl files, %d look like count matrices",
                accession, len(urls), len(count_urls))

    for url in other_urls:
        record["skipped"].append(Path(url).name)

    if not count_urls:
        record["status"] = "no_count_file"
        logger.info("  %s: no count matrix file identified among: %s",
                    accession, [Path(u).name for u in urls[:5]])
        return record

    for url in count_urls:
        name = Path(url).name
        dest = cache_dir / accession / name

        size = file_size(url)
        if size and size > MAX_BYTES:
            logger.info("  Skipping %s (%.0f MB > limit)", name, size / 1e6)
            record["skipped"].append(name)
            continue

        if dry_run:
            logger.info("  [DRY RUN] Would download: %s", name)
            record["downloaded"].append(name)
            continue

        success = download_file(url, dest)
        if success:
            record["downloaded"].append(name)

    record["status"] = "downloaded" if record["downloaded"] else "no_download"
    return record


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accessions", nargs="+",
                        help="Specific GEO accessions to process (default: all candidates)")
    parser.add_argument("--per-study-dir", type=Path,
                        default=REPO_ROOT / "results" / "per_study" / "transcriptomics")
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be downloaded without downloading")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.accessions:
        candidates = args.accessions
    else:
        candidates = find_rna_seq_candidates(args.cache_dir, args.per_study_dir)

    logger.info("Studies to process: %d", len(candidates))

    results = []
    for acc in candidates:
        rec = process_study(acc, args.cache_dir, dry_run=args.dry_run)
        results.append(rec)

    # Summary
    by_status: dict[str, list[str]] = {}
    for r in results:
        by_status.setdefault(r["status"], []).append(r["accession"])

    print("\n========== Download Summary ==========")
    for status, accs in sorted(by_status.items()):
        print(f"  {status}: {len(accs)}")
        if status in ("downloaded", "no_count_file") and len(accs) <= 15:
            for acc in accs:
                r = next(x for x in results if x["accession"] == acc)
                extra = r.get("downloaded") or r.get("skipped", [])[:2]
                print(f"    {acc}: {extra}")


if __name__ == "__main__":
    main()
