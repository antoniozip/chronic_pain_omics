"""Extract and merge per-sample count files from GEO _RAW.tar tarballs (Step 2.1).

For each eligible study that has a groups.csv but no per-study effects file:
  1. Downloads {acc}_RAW.tar from GEO FTP (in-memory; tarballs are <30 MB for bulk RNA-seq)
  2. Skips single-cell studies (barcodes.tsv.gz / *.h5 / matrix.mtx.gz content)
  3. Parses per-sample count files from the tarball:
       - 2-column HTSeq  : gene_id<TAB>count  (skips __xxx summary lines)
       - featureCounts   : GeneID<TAB>count   (extra empty tab columns stripped)
       - RSEM            : uses expected_count column
       - UTF-16 encoded  : auto-detected and decoded
       - Float values    : passed through (limma-trend in downstream R)
  4. Filters to case/control samples defined in {acc}_groups.csv
  5. Names columns with vocabulary tokens (sham_control_N, cfa_case_N, etc.)
  6. Merges all samples into one matrix; saves as {acc}_raw_count_merged.txt.gz

Output format is priority-1 in find_suppl_count_file() (raw_count prefix pattern).

Usage:
  python pipeline/05d_extract_raw_tar.py
  python pipeline/05d_extract_raw_tar.py --accessions GSE272517 GSE138024
  python pipeline/05d_extract_raw_tar.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import os
import re
import tarfile
import tempfile
import time
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)

REPO_ROOT  = Path(__file__).resolve().parent.parent
CACHE_DIR  = REPO_ROOT / "data" / "raw" / "geo_cache"
PER_STUDY  = REPO_ROOT / "results" / "per_study" / "transcriptomics"

GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
MEM_MAX_BYTES  = 100 * 1024 * 1024   # 100 MB — download in-memory below this
DISK_MAX_BYTES = 2 * 1024 * 1024 * 1024  # 2 GB — download to temp file up to this
RATE_DELAY   = 1.0

# ---------------------------------------------------------------------------
# Vocabulary maps: condition string fragment → canonical token
# ---------------------------------------------------------------------------

PAIN_VOCAB: list[tuple[str, str]] = [
    ("paclitaxel", "paclitaxel"),
    ("oxaliplatin", "oxaliplatin"),
    ("cisplatin",   "cisplatin"),
    ("vincristine", "vcr"),
    ("snl",         "sni"),
    ("psnl",        "sni"),
    ("sni",         "sni"),
    ("cci",         "cci"),
    ("cfa",         "cfa"),
    ("mia",         "oa"),
    ("flit",        "injured"),
    ("spared",      "sni"),
    ("ligat",       "ligated"),
    ("inj",         "injured"),
    ("pain",        "pain"),
]

CTRL_VOCAB: list[tuple[str, str]] = [
    ("sham",    "sham"),
    ("naive",   "naive"),
    ("vehicle", "vehicle"),
    ("saline",  "saline"),
    ("normal",  "naive"),
    ("control", "naive"),
    ("ctrl",    "naive"),
    ("veh",     "vehicle"),
    ("con",     "naive"),
]

# Filename fragments that identify non-count supplementary files to skip
SKIP_NAMES = re.compile(
    r"\.(bgx|gtf|gff|bed|bam|bai|fa|fasta|gff3)(\.gz)?$"
    r"|_annotation\.|_genome\.",
    re.IGNORECASE,
)

# Single-cell file signatures — if any member matches, skip the whole study
SCRNA_SIGNATURES = re.compile(
    r"barcodes\.tsv(\.gz)?$|features\.tsv(\.gz)?$"
    r"|matrix\.mtx(\.gz)?$|filtered_feature_bc_matrix",
    re.IGNORECASE,
)

# HTSeq summary lines to skip
HTSEQ_SKIP = re.compile(r"^__[a-z_]+$", re.IGNORECASE)

# Ensembl version suffix e.g. ENSMUSG00000000001.4 → ENSMUSG00000000001
ENSEMBL_VERSION = re.compile(r"^(ENS[A-Z]+G\d+)\.\d+$")


# ---------------------------------------------------------------------------
# URL helpers
# ---------------------------------------------------------------------------

def raw_tar_url(accession: str) -> str:
    n      = accession[3:]           # strip "GSE"
    prefix = n[:-3] + "nnn"
    return f"{GEO_FTP_BASE}/GSE{prefix}/{accession}/suppl/{accession}_RAW.tar"


def fetch_tarball(url: str, max_bytes: int = DISK_MAX_BYTES) -> str | bytes | None:
    """Download a tarball; returns bytes for small ones, tempfile path for large ones, or None."""
    try:
        head = requests.head(url, timeout=15, allow_redirects=True)
        size = head.headers.get("Content-Length")
        if size and int(size) > max_bytes:
            logger.warning("  Tarball too large (%.0f GB > %.0f GB limit) — skipping",
                           int(size) / 1e9, max_bytes / 1e9)
            return None

        if size and int(size) <= MEM_MAX_BYTES:
            # Small enough for in-memory
            resp = requests.get(url, timeout=300)
            resp.raise_for_status()
            logger.info("  Downloaded %.1f MB (in-memory)", len(resp.content) / 1e6)
            return resp.content
        else:
            # Stream to temp file for larger downloads
            suffix = ".tar"
            fd, tmp_path = tempfile.mkstemp(suffix=suffix, prefix="geo_tar_")
            os.close(fd)
            logger.info("  Streaming to %s ...", Path(tmp_path).name)
            with requests.get(url, stream=True, timeout=600) as resp:
                resp.raise_for_status()
                downloaded = 0
                with open(tmp_path, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=65536):
                        fh.write(chunk)
                        downloaded += len(chunk)
            logger.info("  Downloaded %.1f MB (to disk)", downloaded / 1e6)
            return tmp_path
    except Exception as exc:
        logger.warning("  Fetch failed: %s", exc)
        return None


# ---------------------------------------------------------------------------
# File format detection and parsing
# ---------------------------------------------------------------------------

def _decode_bytes(raw: bytes) -> str:
    """Decode bytes, handling UTF-16 LE/BE with BOM."""
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16", errors="replace")
    if raw[:3] == b"\xef\xbb\xbf":
        return raw[3:].decode("utf-8", errors="replace")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1", errors="replace")


def _read_member(tf: tarfile.TarFile, member: tarfile.TarInfo) -> str | None:
    """Extract and decode a tarfile member to string."""
    fobj = tf.extractfile(member)
    if fobj is None:
        return None
    raw = fobj.read()
    # Try gzip decompression first
    if member.name.endswith(".gz"):
        try:
            raw = gzip.decompress(raw)
        except Exception:
            pass
    return _decode_bytes(raw)


def _clean_gene_id(gid: str) -> str:
    """Strip Ensembl version suffix for cleaner IDs."""
    m = ENSEMBL_VERSION.match(gid)
    return m.group(1) if m else gid


def parse_sample_content(text: str) -> pd.Series | None:
    """
    Parse a per-sample count file into a pd.Series (index=gene_id, values=count).

    Handles:
      - HTSeq 2-col (gene_id TAB count, skip __xxx lines)
      - featureCounts (GeneID TAB count TAB TAB ...)
      - RSEM (multi-column with expected_count)
      - plain 2-col (gene TAB value) with optional header
    """
    lines = [line.rstrip("\r\n") for line in text.split("\n") if line.strip()]
    if not lines:
        return None

    # Skip featureCounts/SAM-style comment lines starting with #
    data_lines = [line for line in lines if not line.startswith("#")]
    if not data_lines:
        return None

    first = data_lines[0]
    fields = first.split("\t")

    # --- RSEM detection: header contains "expected_count" ---
    if any("expected_count" in f.lower() for f in fields):
        header = [f.lower() for f in fields]
        ec_idx = next(i for i, h in enumerate(header) if "expected_count" in h)
        # gene ID is first column
        rows: dict[str, float] = {}
        for line in data_lines[1:]:
            parts = line.split("\t")
            if len(parts) <= max(0, ec_idx):
                continue
            gid = _clean_gene_id(parts[0].strip())
            if not gid:
                continue
            try:
                rows[gid] = float(parts[ec_idx])
            except ValueError:
                continue
        return pd.Series(rows, dtype=float) if rows else None

    # --- featureCounts detection: header starts with "Geneid" or "GeneID",
    #     and data lines have trailing empty tab fields ---
    if fields[0].lower() in ("geneid", "gene_id") and len(fields) > 2:
        # Look for a non-empty second field in data
        rows = {}
        for line in data_lines[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gid = _clean_gene_id(parts[0].strip())
            val_str = parts[1].strip()
            if not gid or not val_str:
                continue
            try:
                rows[gid] = float(val_str)
            except ValueError:
                continue
        return pd.Series(rows, dtype=float) if rows else None

    # --- Detect header and find first numeric value column ---
    # Some formats (e.g. Cufflinks FPKM per-gene) have extra annotation columns
    # before the numeric value (Gene_Name, Locus, etc.).
    has_header = False
    if len(fields) >= 2:
        try:
            float(fields[1])
        except ValueError:
            has_header = True

    data_start = 1 if has_header else 0
    # Determine value column index by scanning first few data lines
    value_col = 1
    for probe in data_lines[data_start: data_start + 5]:
        parts_p = probe.split("\t")
        for ci in range(1, len(parts_p)):
            try:
                float(parts_p[ci].strip())
                value_col = ci
                break
            except ValueError:
                continue
        else:
            continue
        break

    rows = {}
    for line in data_lines[data_start:]:
        parts = line.split("\t")
        if len(parts) <= value_col:
            continue
        gid = parts[0].strip()
        if HTSEQ_SKIP.match(gid):
            continue
        gid = _clean_gene_id(gid)
        val_str = parts[value_col].strip()
        if not gid or not val_str:
            continue
        try:
            rows[gid] = float(val_str)
        except ValueError:
            continue
    return pd.Series(rows, dtype=float) if rows else None


# ---------------------------------------------------------------------------
# Sample identification helpers
# ---------------------------------------------------------------------------

DERIVED_QUANT = re.compile(r"tpm|fpkm|rpkm|\bcpm\b|normali[sz]ed", re.IGNORECASE)
RAW_QUANT = re.compile(
    r"raw[._-]?count|htseq|feature[._-]?count|expected[._-]?count|counts?\b",
    re.IGNORECASE,
)


def quant_rank(name: str) -> int:
    """Rank a tar member by how directly it carries counts. Lower is better.

    A name saying both ("..._FPKM_counts.txt.gz") is read as derived: the
    derived label is the one that changes what the numbers mean.
    """
    stem = Path(name).name
    if DERIVED_QUANT.search(stem):
        return 2
    if RAW_QUANT.search(stem):
        return 0
    return 1


def select_one_member_per_sample(
    members: list[tarfile.TarInfo],
) -> list[tarfile.TarInfo]:
    """Keep at most one file per GSM, preferring raw counts.

    GSE205494 deposits an FPKM file *and* a raw-count file per sample. Adding a
    column per file doubled its sample count and mixed two quantification
    scales in one matrix, which `is_integer_counts` then reads as continuous
    and hands to limma-trend as though it were one assay. Members with no GSM
    in their name are passed through untouched, because the auto-detection
    fallback keys on filenames instead.
    """
    best: dict[str, tarfile.TarInfo] = {}
    out: list[tarfile.TarInfo] = []
    for member in members:
        gsm = gsm_from_filename(member.name)
        if gsm is None:
            out.append(member)
            continue
        held = best.get(gsm)
        if held is None or (quant_rank(member.name), member.name) < (
            quant_rank(held.name),
            held.name,
        ):
            best[gsm] = member
    out.extend(best[g] for g in sorted(best))
    return out


def gsm_from_filename(name: str) -> str | None:
    """Extract GSM accession from a filename like GSM12345_title.txt.gz."""
    m = re.match(r"^(GSM\d+)", Path(name).name, re.IGNORECASE)
    return m.group(1) if m else None


def condition_from_filename(name: str) -> str:
    """Extract the condition part from GSM12345_condition.ext.gz."""
    stem = Path(name).name
    # Strip .gz and remaining extension
    for _ in range(3):
        stem = Path(stem).stem if "." in stem else stem
    # Strip leading GSM accession
    stem = re.sub(r"^GSM\d+[_\-]?", "", stem, flags=re.IGNORECASE)
    return stem


def make_col_name(condition: str, group: str, idx: int) -> str:
    """Create a vocabulary-matchable column name."""
    cond_lo = condition.lower().replace("-", "_").replace(" ", "_")

    if group == "case":
        for frag, token in PAIN_VOCAB:
            if frag in cond_lo:
                return f"{token}_case_{idx}"
        return f"injured_case_{idx}"
    else:
        for frag, token in CTRL_VOCAB:
            if frag in cond_lo:
                return f"{token}_control_{idx}"
        return f"naive_control_{idx}"


# ---------------------------------------------------------------------------
# Per-study extraction
# ---------------------------------------------------------------------------

def load_groups(accession: str, cache_dir: Path) -> dict[str, str]:
    """Return {gsm_id: group} from {accession}_groups.csv."""
    path = cache_dir / f"{accession}_groups.csv"
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {r["sample_id"]: r["group"] for r in csv.DictReader(fh)}


def is_sc_study(members: list[tarfile.TarInfo]) -> bool:
    """Return True if tarball appears to contain single-cell data."""
    for m in members:
        if SCRNA_SIGNATURES.search(m.name):
            return True
    return False


def _discard_tarball(tmp_path: str | None) -> None:
    """Remove a streamed tarball once nothing else will re-open it."""
    if tmp_path and os.path.exists(tmp_path):
        os.unlink(tmp_path)


def extract_study(accession: str, cache_dir: Path, dry_run: bool = False) -> str:
    """
    Extract, merge, and save the count matrix for one study.
    Returns a status string.
    """
    out_dir = cache_dir / accession
    out_path = out_dir / f"{accession}_raw_count_merged.txt.gz"
    if out_path.exists():
        logger.info("%s: merged matrix already present — skipping", accession)
        return "already_present"

    groups = load_groups(accession, cache_dir)
    if not groups:
        logger.warning("%s: no groups.csv found — skipping", accession)
        return "no_groups"

    url = raw_tar_url(accession)
    logger.info("%s: downloading %s", accession, url.split("/")[-1])

    tarball = fetch_tarball(url)
    if tarball is None:
        return "download_failed"

    tmp_path: str | None = None
    try:
        if isinstance(tarball, bytes):
            tf = tarfile.open(fileobj=io.BytesIO(tarball))
        else:
            tmp_path = tarball  # str path to temp file
            tf = tarfile.open(name=tmp_path)
    except Exception as exc:
        logger.warning("%s: cannot open tarball: %s", accession, exc)
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)
        return "tar_error"

    members = tf.getmembers()
    logger.info("%s: tarball contains %d files", accession, len(members))

    if is_sc_study(members):
        logger.info("%s: single-cell data detected — skipping", accession)
        tf.close()
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)
        return "scrna_skipped"

    # Filter to count-like members only
    count_members = [
        m for m in members
        if m.isfile()
        and not SKIP_NAMES.search(m.name)
        and m.name.lower().endswith((".txt.gz", ".txt", ".csv.gz", ".csv",
                                     ".tsv.gz", ".tsv"))
    ]
    # One file per GSM. A tarball carrying both an FPKM and a raw-count file
    # per sample would otherwise contribute two columns for one animal.
    n_before = len(count_members)
    count_members = select_one_member_per_sample(count_members)
    if len(count_members) < n_before:
        logger.info("%s: %d of %d files dropped as duplicate quantifications",
                    accession, n_before - len(count_members), n_before)
    logger.info("%s: %d candidate count files", accession, len(count_members))

    case_idx   = 0
    ctrl_idx   = 0
    series_map: dict[str, pd.Series] = {}

    for member in count_members:
        gsm = gsm_from_filename(member.name)
        if gsm is None:
            logger.debug("  Skipping (no GSM ID): %s", member.name)
            continue

        group = groups.get(gsm, "")
        if group not in ("case", "control"):
            logger.debug("  %s (%s): group=%r — excluded", member.name, gsm, group)
            continue

        text = _read_member(tf, member)
        if text is None:
            logger.warning("  %s: could not read", member.name)
            continue

        series = parse_sample_content(text)
        if series is None or series.empty:
            logger.warning("  %s: parsing returned empty series", member.name)
            continue

        condition = condition_from_filename(member.name)
        if group == "case":
            case_idx += 1
            col_name = make_col_name(condition, "case", case_idx)
        else:
            ctrl_idx += 1
            col_name = make_col_name(condition, "control", ctrl_idx)

        logger.info("  %s → %s (%d genes)", gsm, col_name, len(series))
        series_map[col_name] = series

    tf.close()
    # The temp file is NOT removed here. Both auto-detection fallbacks below
    # re-open the tarball, and deleting it first made every study that needed
    # them die on its own download with FileNotFoundError -- GSE255553 and
    # GSE289659 in the 2026-08-29 amendment run. Cleanup happens at each exit
    # instead, via _discard_tarball().

    n_case = sum(1 for k in series_map if "_case_" in k)
    n_ctrl = sum(1 for k in series_map if "_control_" in k)
    logger.info("%s: %d case, %d control columns parsed (groups-based)", accession, n_case, n_ctrl)

    # Fallback: if groups didn't match any tarball files, try auto-detection
    if n_case == 0 and n_ctrl == 0:
        logger.info("%s: groups.csv GSM IDs didn't match tarball — trying auto-detection",
                    accession)
        # Re-open tarball for fallback pass
        if isinstance(tarball, bytes):
            tf2 = tarfile.open(fileobj=io.BytesIO(tarball))
        elif tmp_path and os.path.exists(tmp_path):
            tf2 = tarfile.open(name=tmp_path)
        else:
            tf2 = tarfile.open(name=tarball)  # type: ignore[arg-type]
        members2 = tf2.getmembers()
        count_members2 = [
            m for m in members2
            if m.isfile()
            and not SKIP_NAMES.search(m.name)
            and m.name.lower().endswith((".txt.gz", ".txt", ".csv.gz", ".csv",
                                         ".tsv.gz", ".tsv"))
        ]
        count_members2 = select_one_member_per_sample(count_members2)
        case_idx = ctrl_idx = 0
        series_map.clear()

        for member in count_members2:
            condition = condition_from_filename(member.name)
            cond_lo = condition.lower().replace("-", "_").replace(" ", "_")

            # Check pain vocab
            is_case = any(frag in cond_lo for frag, _ in PAIN_VOCAB)
            is_ctrl = any(frag in cond_lo for frag, _ in CTRL_VOCAB)

            if is_case and not is_ctrl:
                group = "case"
            elif is_ctrl and not is_case:
                group = "control"
            elif is_case and is_ctrl:
                continue  # ambiguous
            else:
                continue  # can't determine

            text = _read_member(tf2, member)
            if text is None:
                continue
            series = parse_sample_content(text)
            if series is None or series.empty:
                continue

            if group == "case":
                case_idx += 1
                col_name = make_col_name(condition, "case", case_idx)
            else:
                ctrl_idx += 1
                col_name = make_col_name(condition, "control", ctrl_idx)

            logger.info("  [auto] %s → %s (%d genes)",
                        member.name.split("/")[-1], col_name, len(series))
            series_map[col_name] = series

        # Second pass: if still nothing or too few, try by-exclusion strategy
        if n_case := (sum(1 for k in series_map if "_case_" in k) < 2
                      or sum(1 for k in series_map if "_control_" in k) < 2):
            logger.info("%s: vocab-only auto-detection insufficient — trying by-exclusion",
                        accession)
            # Re-classify: any non-control file is case
            case_idx2 = ctrl_idx2 = 0
            series_map.clear()
            tf3 = (tarfile.open(fileobj=io.BytesIO(tarball))
                   if isinstance(tarball, bytes)
                   else tarfile.open(
                       name=tarball if isinstance(tarball, str) else tmp_path or ""))
            count_members3 = [
                m for m in tf3.getmembers()
                if m.isfile()
                and not SKIP_NAMES.search(m.name)
                and m.name.lower().endswith((".txt.gz", ".txt", ".csv.gz", ".csv",
                                             ".tsv.gz", ".tsv"))
            ]
            count_members3 = select_one_member_per_sample(count_members3)
            for member in count_members3:
                condition = condition_from_filename(member.name)
                cond_lo = condition.lower().replace("-", "_").replace(" ", "_")
                is_ctrl = any(frag in cond_lo for frag, _ in CTRL_VOCAB)
                group = "control" if is_ctrl else "case"

                text = _read_member(tf3, member)
                if text is None:
                    continue
                series = parse_sample_content(text)
                if series is None or series.empty:
                    continue

                if group == "case":
                    case_idx2 += 1
                    col_name = make_col_name(condition, "case", case_idx2)
                else:
                    ctrl_idx2 += 1
                    col_name = make_col_name(condition, "control", ctrl_idx2)
                series_map[col_name] = series
            tf3.close()

        tf2.close()

    n_case = sum(1 for k in series_map if "_case_" in k)
    n_ctrl = sum(1 for k in series_map if "_control_" in k)
    logger.info("%s: %d case, %d control columns parsed", accession, n_case, n_ctrl)

    if n_case < 2 or n_ctrl < 2:
        logger.warning("%s: insufficient samples (need ≥2/group) — not saving", accession)
        _discard_tarball(tmp_path)
        return "insufficient_samples"

    if dry_run:
        logger.info("  [DRY RUN] Would save %s", out_path)
        _discard_tarball(tmp_path)
        return "dry_run"

    # Merge: outer join — missing values become NaN (R handles them)
    mat = pd.concat(series_map, axis=1)
    mat.index.name = "gene_id"
    mat = mat.reset_index()

    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8") as fh:
        mat.to_csv(fh, sep="\t", index=False)

    logger.info("%s: saved %s (%d genes × %d samples)",
                accession, out_path.name, len(mat), len(series_map))
    _discard_tarball(tmp_path)
    return "saved"


# ---------------------------------------------------------------------------
# Study selection
# ---------------------------------------------------------------------------

def find_raw_tar_candidates(cache_dir: Path, per_study_dir: Path) -> list[str]:
    """Return studies with groups.csv but no effects file and a _RAW.tar on GEO."""
    annotated = {p.stem.replace("_groups", "") for p in cache_dir.glob("*_groups.csv")}
    processed = set()
    for p in per_study_dir.glob("*_effects*.csv"):
        if "normalized_all" not in p.name:
            processed.add(p.stem.split("_effects")[0])

    # Also skip studies that already have any count file in their subdir
    skipped: set[str] = set()
    for acc in annotated - processed:
        sub = cache_dir / acc
        if sub.exists() and any(
            re.search(r"(raw.count|count|matrix|tpm|fpkm)", f.name, re.IGNORECASE)
            for f in sub.iterdir()
            if f.is_file()
        ):
            skipped.add(acc)

    candidates = sorted((annotated - processed) - skipped)

    # Filter to studies that actually have a _RAW.tar (check via HEAD)
    tar_candidates = []
    for acc in candidates:
        url = raw_tar_url(acc)
        try:
            r = requests.head(url, timeout=10, allow_redirects=True)
            if r.status_code == 200:
                tar_candidates.append(acc)
        except Exception:
            pass
        time.sleep(0.3)

    return tar_candidates


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accessions", nargs="+",
                        help="Specific accessions (default: auto-discover)")
    parser.add_argument("--cache-dir",  type=Path, default=CACHE_DIR)
    parser.add_argument("--per-study-dir", type=Path, default=PER_STUDY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.accessions:
        candidates = args.accessions
    else:
        logger.info("Discovering RAW.tar candidates (HEAD requests) ...")
        candidates = find_raw_tar_candidates(args.cache_dir, args.per_study_dir)

    logger.info("Processing %d studies: %s", len(candidates), candidates)

    results: dict[str, str] = {}
    for acc in candidates:
        logger.info("\n=== %s ===", acc)
        status = extract_study(acc, args.cache_dir, dry_run=args.dry_run)
        results[acc] = status
        time.sleep(RATE_DELAY)

    print("\n========== Extraction Summary ==========")
    by_status: dict[str, list[str]] = {}
    for acc, st in results.items():
        by_status.setdefault(st, []).append(acc)
    for st, accs in sorted(by_status.items()):
        print(f"  {st}: {len(accs)}")
        if len(accs) <= 20:
            for a in accs:
                print(f"    {a}")


if __name__ == "__main__":
    main()
