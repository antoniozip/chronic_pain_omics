#!/usr/bin/env python3
"""
Proteomics Differential Abundance Analysis Pipeline
===================================================
Part of the CP_multi-omics chronic pain meta-analysis.

This script handles the complete workflow:
1. Download quantification files from PRIDE for pain-relevant studies
2. Parse MaxQuant output (proteinGroups.txt) or other quant formats
3. Build protein abundance matrices
4. Perform differential abundance analysis (limma, t-test, etc.)
5. Export results for meta-analysis

Usage:
    .venv/bin/python scripts/proteomics_da_pipeline.py \
        [--accession PXD013362] [--download] [--analyze]

Requires:
    - numpy, pandas, scipy, statsmodels
    - For MaxQuant parsing: pyopenms (optional), pandas
    - For raw reprocessing: MaxQuant (commercial) or FragPipe/MSFragger (free)
"""

from __future__ import annotations

import argparse
import json
import logging
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np

try:
    import pandas as pd
except ImportError:
    pd = None
    print("WARNING: pandas not installed. Install with: pip install pandas openpyxl")

try:
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
except ImportError:
    stats = None
    print("WARNING: scipy/statsmodels not installed. Install with: pip install scipy statsmodels")


# ─── Paths ───────────────────────────────────────────────────────────────────
WORKDIR = Path("/media/antonio/data/CP_multi-omics")
RAW_DIR = WORKDIR / "data/raw/proteomics"
INTERIM_DIR = WORKDIR / "data/interim/proteomics"
RESULTS_DIR = WORKDIR / "results/proteomics"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PRIDE_FTP_BASE = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive"

logger = logging.getLogger("proteomics_da")


# ─── Configuration ───────────────────────────────────────────────────────────
# Pain-relevant datasets identified from metadata screening
PAIN_RELEVANT_DATASETS = [
    {
        "accession": "PXD013362",
        "title": "PACAP and other neuropeptide targets link chronic migraine "
                 "and opioid-induced hyperalgesia in mouse models",
        "species": "Mus musculus (mouse)",
        "relevance_score": 40,
        "matched_terms": ["PACAP", "chronic migraine", "hyperalgesia", "migraine", "opioid"],
        "quant_files": [
            "peptide_quantitation.xlsx",       # ~1 MB, easiest
            "Quantitation_cohort_1.rar",       # 315 GB, large archive
            "Quantitation_cohort_2.rar",       # 452 GB, large archive
        ],
        "raw_files_count": 14,
        "has_quant_xlsx": True,
        "pubmed_id": None,  # TODO: look up from PRIDE metadata
        "groups": ["wildtype", "PACAP_KO"],    # inferred from title
        "tissue": ["hypothalamus", "PAG", "NAc", "RVM", "dorsal_horn", "trigeminal"],
    },
    {
        "accession": "PXD015949",
        "title": "Quantitative characterization of the neuropeptide level "
                 "changes in dorsal horn and dorsal root ganglia regions of "
                 "the murine itch models",
        "species": "Mus musculus (mouse)",
        "relevance_score": 8,
        "matched_terms": ["dorsal horn", "dorsal root"],
        "quant_files": [
            "peptides_DRG.mzid.gz",
            "peptides_1_1_0.mzid.gz",
        ],
        "raw_files_count": 3,
        "has_quant_xlsx": False,
        "pubmed_id": None,
        "groups": ["control", "itch_model"],
        "tissue": ["dorsal_horn", "dorsal_root_ganglion"],
    },
]

# DA analysis parameters
DA_CONFIG = {
    "padj_threshold": 0.05,
    "lfc_threshold": 1.0,
    "min_samples_per_group": 2,
    "min_peptides_per_protein": 1,
    "imputation_method": "min_detected",  # or "knn", "bpca", "none"
    "normalization": "LFQ",              # LFQ intensity, TMT, or spectral counts
}


# ─── Utility Functions ──────────────────────────────────────────────────────

def create_ssl_context() -> ssl.SSLContext:
    """Create a permissive SSL context for FTP/HTTPS access to public data."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def pride_ftp_url(accession: str, filename: str) -> str:
    """Construct PRIDE FTP URL for a file."""
    # Submission date determines path; default to pattern
    return f"{PRIDE_FTP_BASE}/2019/11/{accession}/{filename}"


def pride_file_list_url(accession: str) -> str:
    return f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files"


def fetch_pride_file_list(accession: str) -> list[dict]:
    """Get file listing including FTP paths for a PRIDE project."""
    ctx = create_ssl_context()
    url = pride_file_list_url(accession)
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, context=ctx, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        file_list = data if isinstance(data, list) else \
                    data.get("_embedded", {}).get("files", []) if isinstance(data, dict) else []
        return file_list
    except Exception as e:
        logger.error(f"Failed to fetch file list for {accession}: {e}")
        return []


def get_ftp_urls_for_file(file_entry: dict) -> list[str]:
    """Extract FTP URLs from a PRIDE file entry."""
    urls = []
    locations = file_entry.get("publicFileLocations", [])
    for loc in locations:
        if isinstance(loc, dict):
            name = loc.get("name", "")
            value = loc.get("value", "")
            if "FTP" in name:
                urls.append(value)
    return urls


# ─── Step 1: Download Quant Files ──────────────────────────────────────────

def download_file(url: str, dest: Path, max_retries: int = 3, chunk_size: int = 8192) -> bool:
    """Download a file with retries, showing progress for large files."""
    logger.info(f"Downloading {url} -> {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=300) as resp:
                total = int(resp.headers.get("content-length", 0))
                downloaded = 0
                with open(dest, "wb") as f:
                    while chunk := resp.read(chunk_size):
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total and downloaded % (10 * 1024 * 1024) == 0:
                            pct = downloaded / total * 100
                            logger.info(
                                f"  {downloaded//1024//1024}MB / "
                                f"{total//1024//1024}MB ({pct:.0f}%)")
            logger.info(f"  Downloaded {dest} ({dest.stat().st_size / 1024 / 1024:.1f} MB)")
            return True
        except Exception as e:
            logger.warning(f"  Attempt {attempt+1} failed: {e}")
            if attempt < max_retries - 1:
                time.sleep(5)
    logger.error(f"  Failed to download {url} after {max_retries} attempts")
    return False


def download_all_quant_files(accession: str, file_list: list[dict], output_dir: Path) -> list[Path]:
    """Download quantification-related files for a PRIDE project."""
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    quant_extensions = [".xlsx", ".csv", ".tsv", ".txt", ".mzid.gz", ".mztab.gz"]
    quant_keywords = ["quant", "protein", "peptide", "intensity", "maxquant", "proteingroups"]

    for f_entry in file_list:
        fname = f_entry.get("fileName", "")
        fname_lower = fname.lower()

        # Check if this is a quant file
        is_quant = any(fname_lower.endswith(ext) for ext in quant_extensions)
        is_quant = is_quant or any(kw in fname_lower for kw in quant_keywords)

        # Skip raw files and MGF
        if fname_lower.endswith((".raw", ".mgf", ".wiff", ".mzml", ".mzxml")):
            continue

        if not is_quant:
            continue

        # Get FTP URL
        ftp_urls = get_ftp_urls_for_file(f_entry)
        if not ftp_urls:
            continue

        dest = output_dir / fname
        if dest.exists() and dest.stat().st_size > 0:
            logger.info(f"  Already exists: {dest} ({dest.stat().st_size / 1024:.1f} KB)")
            downloaded.append(dest)
            continue

        # Try each URL
        for url in ftp_urls:
            if download_file(url, dest):
                downloaded.append(dest)
                break

    return downloaded


# ─── Step 2: Parse Quantification Data ──────────────────────────────────────

def parse_protein_groups_txt(path: Path) -> pd.DataFrame | None:
    """Parse a MaxQuant proteinGroups.txt file."""
    if pd is None:
        logger.error("pandas required for parsing")
        return None
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        logger.info(f"  Parsed {path}: {df.shape[0]} proteins, {df.shape[1]} columns")
        return df
    except Exception as e:
        logger.error(f"  Failed to parse {path}: {e}")
        return None


def parse_quant_xlsx(path: Path) -> pd.DataFrame | None:
    """Parse an Excel quantification file."""
    if pd is None:
        logger.error("pandas required for parsing")
        return None
    try:
        xl = pd.ExcelFile(path)
        logger.info(f"  Excel sheets: {xl.sheet_names}")
        # Try to find the main data sheet
        for sheet in xl.sheet_names:
            if any(kw in sheet.lower() for kw in
                   ["protein", "peptide", "quant", "data", "intensity"]):
                df = pd.read_excel(path, sheet_name=sheet)
                logger.info(f"  Parsed sheet '{sheet}': {df.shape}")
                return df
        # Fall back to first sheet
        df = pd.read_excel(path, sheet_name=xl.sheet_names[0])
        logger.info(f"  Parsed first sheet '{xl.sheet_names[0]}': {df.shape}")
        return df
    except Exception as e:
        logger.error(f"  Failed to parse Excel {path}: {e}")
        return None


def parse_mzid(path: Path) -> pd.DataFrame | None:
    """Parse an mzIdentML file - requires pyopenms or similar."""
    # This is a placeholder - mzIdentML parsing needs specialized tools
    logger.info(f"  mzIdentML parsing not implemented yet for {path.name}")
    logger.info("  Recommend using ProteoWizard or pyopenms for mzid parsing")
    return None


def parse_quant_file(path: Path) -> pd.DataFrame | None:
    """Auto-detect format and parse a quantification file."""
    ext = path.suffix.lower()
    if ext == ".xlsx":
        return parse_quant_xlsx(path)
    elif path.name.endswith("proteinGroups.txt"):
        return parse_protein_groups_txt(path)
    elif path.name.endswith(".txt") or ext == ".tsv":
        if pd is not None:
            try:
                df = pd.read_csv(path, sep="\t", low_memory=False)
                return df
            except (OSError, UnicodeDecodeError, pd.errors.ParserError,
                    pd.errors.EmptyDataError) as exc:
                logger.warning("could not parse %s as TSV: %s", path.name, exc)
    elif path.name.endswith(".csv"):
        if pd is not None:
            try:
                df = pd.read_csv(path, low_memory=False)
                return df
            except (OSError, UnicodeDecodeError, pd.errors.ParserError,
                    pd.errors.EmptyDataError) as exc:
                logger.warning("could not parse %s as CSV: %s", path.name, exc)
    elif path.name.endswith((".mzid.gz", ".mzid")):
        return parse_mzid(path)
    return None


# ─── Step 3: Build Protein Abundance Matrix ─────────────────────────────────

def build_abundance_matrix(
    parsed_data: pd.DataFrame,
    accession: str,
    intensity_columns: list[str] | None = None,
    protein_id_col: str | None = None,
) -> pd.DataFrame:
    """
    Build a protein abundance matrix suitable for differential analysis.
    
    Expected format from MaxQuant: columns like 'LFQ intensity Sample1', 'LFQ intensity Sample2'
    Or from generic quant: columns like 'Sample1', 'Sample2' with protein IDs as index.
    """
    if parsed_data is None:
        return pd.DataFrame()

    # Try to identify columns
    if intensity_columns is None:
        # Auto-detect: look for LFQ intensity, iBAQ, or raw intensity columns
        intensity_cols = [c for c in parsed_data.columns 
                          if any(x in c.lower() for x in
                                 ["lfq intensity", "lfq", "ibaq",
                                  "intensity "])]
        if not intensity_cols:
            # Fall back to columns that look like sample names (not metadata)
            meta_cols = ["protein ids", "protein id", "gene name", "gene names", "peptides", 
                         "peptide counts", "score", "coverage", "q-value", "sequence coverage"]
            intensity_cols = [c for c in parsed_data.columns 
                              if not any(m in c.lower() for m in meta_cols)]
    else:
        intensity_cols = intensity_columns

    if not intensity_cols:
        logger.warning(f"  No intensity columns found in {accession} data")
        return pd.DataFrame()

    # Try to find protein ID column
    if protein_id_col is None:
        for candidate in ["Protein IDs", "Protein ID", "protein_ids", "protein_id", 
                          "protein", "Majority protein IDs", "Accession"]:
            if candidate in parsed_data.columns:
                protein_id_col = candidate
                break

    if protein_id_col:
        matrix = parsed_data[[protein_id_col] + intensity_cols].copy()
        matrix = matrix.set_index(protein_id_col)
    else:
        matrix = parsed_data[intensity_cols].copy()
        # Use index as protein IDs
        matrix.index = [f"{accession}_protein_{i}" for i in range(len(matrix))]

    # Log2 transform
    matrix = matrix.replace(0, np.nan)
    matrix_log2 = np.log2(matrix)

    logger.info(f"  Built abundance matrix: {matrix_log2.shape} ({len(intensity_cols)} samples)")
    return matrix_log2


# ─── Step 4: Differential Abundance Analysis ────────────────────────────────

def differential_abundance(
    matrix: pd.DataFrame,
    group_labels: list[str] | None = None,
    case_indices: list[int] | None = None,
    control_indices: list[int] | None = None,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Perform differential abundance analysis using Welch's t-test.
    
    For more sophisticated analysis (limma, moderated t-test), the
    matrix should be exported to R or use the statsmodels framework.
    
    Parameters:
    -----------
    matrix: log2-transformed abundance matrix (proteins x samples)
    group_labels: list of case/control per sample column
    case_indices: indices of case samples
    control_indices: indices of control samples
    
    Returns:
    --------
    DataFrame with log2FC, p-value, adjusted p-value per protein
    """
    if matrix.empty:
        return pd.DataFrame()

    n_proteins, n_samples = matrix.shape

    if group_labels:
        unique_groups = list(set(group_labels))
        if len(unique_groups) < 2:
            logger.error("Need at least 2 groups for DA analysis")
            return pd.DataFrame()
        # Assume first unique group is control, last is case
        control_label = unique_groups[0]
        case_label = unique_groups[-1]
        control_indices = [i for i, g in enumerate(group_labels) if g == control_label]
        case_indices = [i for i, g in enumerate(group_labels) if g == case_label]
    elif case_indices is None or control_indices is None:
        logger.error("Must provide either group_labels or case/control indices")
        return pd.DataFrame()

    logger.info(f"  DA analysis: {len(case_indices)} case vs "
                f"{len(control_indices)} control samples")

    results = []
    for protein_name in matrix.index:
        case_values = matrix.loc[protein_name, matrix.columns[case_indices]].dropna().values
        control_values = matrix.loc[protein_name, matrix.columns[control_indices]].dropna().values

        if len(case_values) < 2 or len(control_values) < 2:
            continue

        mean_case = np.mean(case_values)
        mean_control = np.mean(control_values)
        log2fc = mean_case - mean_control

        try:
            t_stat, p_val = stats.ttest_ind(case_values, control_values, equal_var=False)
        except (ValueError, ZeroDivisionError):
            continue

        results.append({
            "protein_id": protein_name,
            "log2FC": log2fc,
            "mean_case": mean_case,
            "mean_control": mean_control,
            "n_case": len(case_values),
            "n_control": len(control_values),
            "t_statistic": t_stat,
            "p_value": p_val,
        })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)

    # Multiple testing correction
    if len(result_df) > 1:
        reject, padj, _, _ = multipletests(result_df["p_value"], method="fdr_bh")
        result_df["padj"] = padj
    else:
        result_df["padj"] = result_df["p_value"]

    # Call significant
    result_df["significant"] = (result_df["padj"] < padj_threshold) & \
                                (result_df["log2FC"].abs() > lfc_threshold)

    result_df = result_df.sort_values("padj")
    logger.info(f"  Found {result_df['significant'].sum()} significant "
                f"proteins (padj<{padj_threshold}, "
                f"|log2FC|>{lfc_threshold})")

    return result_df


# ─── Step 5: Results Export ─────────────────────────────────────────────────

def export_results(
    result_df: pd.DataFrame,
    accession: str,
    output_dir: Path,
    matrix: pd.DataFrame | None = None,
) -> Path:
    """Save DA results and supporting data."""
    study_dir = output_dir / accession
    study_dir.mkdir(parents=True, exist_ok=True)

    # DA results
    result_path = study_dir / "differential_abundance.csv"
    result_df.to_csv(result_path, index=False)
    logger.info(f"  Saved DA results: {result_path}")

    # Significant proteins only
    if "significant" in result_df.columns:
        sig = result_df[result_df["significant"]]
        sig_path = study_dir / "significant_proteins.csv"
        sig.to_csv(sig_path, index=False)
        logger.info(f"  Saved {len(sig)} significant proteins: {sig_path}")

    # Volcano plot data
    volcano = result_df[["protein_id", "log2FC", "padj"]].copy()
    volcano["-log10_padj"] = -np.log10(volcano["padj"].clip(lower=1e-300))
    volcano_path = study_dir / "volcano_data.csv"
    volcano.to_csv(volcano_path, index=False)
    logger.info(f"  Saved volcano plot data: {volcano_path}")

    # Abundance matrix (if provided)
    if matrix is not None:
        matrix_path = study_dir / "abundance_matrix_log2.csv"
        matrix.to_csv(matrix_path)
        logger.info(f"  Saved abundance matrix: {matrix_path}")

    return result_path


# ─── Full Pipeline ──────────────────────────────────────────────────────────

def run_pipeline(
    accession: str,
    download: bool = False,
    analyze: bool = True,
) -> dict[str, Any]:
    """Run the full proteomics DA pipeline for one dataset."""
    logger.info(f"{'='*60}")
    logger.info(f"Processing: {accession}")
    logger.info(f"{'='*60}")

    result = {"accession": accession, "status": "unknown", "steps": {}}

    # Find dataset config
    ds_config = next((d for d in PAIN_RELEVANT_DATASETS if d["accession"] == accession), None)
    if not ds_config:
        logger.warning(f"  No configuration found for {accession}")
        result["status"] = "no_config"
        return result

    study_data_dir = RAW_DIR / accession
    study_data_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Download + discover files
    if download:
        logger.info(f"\n  [Step 1] Downloading files for {accession}")
        file_list = fetch_pride_file_list(accession)
        result["steps"]["file_discovery"] = {"files_found": len(file_list)}

        if not file_list:
            logger.warning("  No files found via API")
            result["status"] = "no_files"
            return result

        downloaded = download_all_quant_files(accession, file_list, study_data_dir)
        result["steps"]["download"] = {"files_downloaded": len(downloaded),
                                        "paths": [str(p) for p in downloaded]}
    else:
        # Just discover what's available
        file_list = fetch_pride_file_list(accession)
        result["steps"]["file_discovery"] = {"files_found": len(file_list)}

        # Check what files already exist locally
        existing_quant = list(study_data_dir.glob("*.xlsx")) + \
                         list(study_data_dir.glob("*.csv")) + \
                         list(study_data_dir.glob("*proteinGroups.txt"))
        downloaded = existing_quant
        if existing_quant:
            logger.info(f"  Found {len(existing_quant)} existing quant files locally")

    # Step 2: Parse quant files
    logger.info("\n  [Step 2] Parsing quantification data")
    parsed_data = None
    for f in downloaded:
        if f.suffix == ".xlsx" and ds_config.get("has_quant_xlsx"):
            parsed_data = parse_quant_xlsx(f)
            if parsed_data is not None:
                logger.info(f"  Successfully parsed {f.name}")
                break
        elif "peptide_quantitation" in f.name:
            parsed_data = parse_quant_xlsx(f)
            if parsed_data is not None:
                break
        elif f.name.endswith("proteinGroups.txt"):
            parsed_data = parse_protein_groups_txt(f)
            if parsed_data is not None:
                break

    if parsed_data is None:
        logger.warning(f"  No parsable quant data found for {accession}")
        result["status"] = "no_parsable_data"
        result["steps"]["parsing"] = {"status": "failed", "note": "No parsable quant files"}
        return result

    result["steps"]["parsing"] = {"format": str(type(parsed_data).__name__),
                                   "shape": list(parsed_data.shape)}

    # Step 3: Build abundance matrix
    logger.info("\n  [Step 3] Building abundance matrix")
    matrix = build_abundance_matrix(parsed_data, accession)
    if matrix.empty:
        logger.warning("  Could not build abundance matrix")
        result["status"] = "matrix_failed"
        return result
    result["steps"]["matrix"] = {"shape": list(matrix.shape)}

    # Step 4: Differential abundance
    if analyze:
        logger.info("\n  [Step 4] Differential abundance analysis")
        # For PXD013362, we know the groups from the study design
        # This needs manual curation - groups are not in the API metadata
        if accession == "PXD013362":
            # The xlsx needs human inspection to determine column-to-group mapping
            logger.info(f"  Group assignment needs manual curation for {accession}")
            logger.info("  See plan/proteomics_analysis_plan.md for details")

        result["status"] = "ready_for_analysis"
    else:
        result["status"] = "data_prepared"

    # Step 5: Export
    if analyze and matrix is not None:
        logger.info("\n  [Step 5] Exporting results")
        # Export intermediate data for manual review
        matrix.to_csv(study_data_dir / "abundance_matrix_log2.csv")
        logger.info("  Saved intermediate abundance matrix to "
                    f"{study_data_dir / 'abundance_matrix_log2.csv'}")

    return result


# ─── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Proteomics differential abundance pipeline for PRIDE datasets"
    )
    parser.add_argument(
        "--accession", "-a",
        default="PXD013362",
        help="PRIDE accession to process (default: PXD013362)"
    )
    parser.add_argument(
        "--download", "-d",
        action="store_true",
        help="Download quantification files from PRIDE"
    )
    parser.add_argument(
        "--analyze", 
        action="store_true",
        help="Run differential abundance analysis"
    )
    parser.add_argument(
        "--list-studies",
        action="store_true",
        help="List available pain-relevant studies"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Process all pain-relevant datasets"
    )
    parser.add_argument(
        "--discover",
        action="store_true",
        help="Discover files for a dataset via PRIDE API without downloading"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose logging"
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.INFO if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.list_studies:
        print("\nPain-Relevant Proteomics Studies from PRIDE:")
        print("=" * 80)
        for ds in PAIN_RELEVANT_DATASETS:
            print(f"  {ds['accession']} (score={ds['relevance_score']})")
            print(f"    {ds['title'][:100]}")
            print(f"    Species: {ds['species']}")
            print(f"    Quant files: {ds['quant_files']}")
            print()
        return

    if args.all:
        results = []
        for ds in PAIN_RELEVANT_DATASETS:
            r = run_pipeline(
                ds["accession"],
                download=args.download,
                analyze=args.analyze,
            )
            results.append(r)
        # Summary
        print("\n\nPipeline Summary:")
        print("=" * 60)
        for r in results:
            print(f"  {r['accession']}: {r['status']}")
        return

    # Single dataset
    result = run_pipeline(
        args.accession,
        download=args.download,
        analyze=args.analyze,
    )

    print(f"\n\nResult for {args.accession}: {result['status']}")
    for step, info in result.get("steps", {}).items():
        print(f"  {step}: {info}")


if __name__ == "__main__":
    main()
