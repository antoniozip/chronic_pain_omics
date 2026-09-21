"""Preprocess count files for studies with non-standard column naming (Step 2.4).

For studies where count file columns don't match GEO sample titles or vocab tokens,
creates renamed/filtered copies with column names that the DA pipeline can recognize.

Studies handled:
  GSE143895 — CCI/SHAM dhsc01_abundance.tsv columns; restore GSM-ID groups.csv
  GSE113941 — CIPN vs Normal input fraction; rename columns to vocab-matchable
  GSE306455 — SNI/LDN/Sham abbreviations (N/L/S); extract SNI vs Sham, rename
  GSE289097 — CFA vs CON (C1/N1); rename to vocab-matchable
  GSE162284 — FLIT vs Sham (Day0/Day3); extract Day3 FLIT vs Day0 Sham, rename
  GSE276193 — 18 sample columns named EM1/CON1 among 14 annotation columns;
              rename to GSM ids using GEO's own per-sample filenames
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

GEO_CACHE = Path("data/raw/geo_cache")


def save_gz(df: pd.DataFrame, path: Path, sep: str = "\t") -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        df.to_csv(fh, sep=sep, index=False)
    logger.info("  Saved %s (%d genes × %d samples)", path.name, len(df), df.shape[1] - 1)


# ---------------------------------------------------------------------------
# GSE143895 — CCI vs SHAM (3+2 samples, DHSC spinal cord)
# ---------------------------------------------------------------------------

def fix_gse143895() -> None:
    logger.info("\n=== GSE143895 (CCI vs SHAM) ===")
    # Restore groups.csv with GSM IDs (auto-annotation correctly assigned these)
    # GSM4276356=dhsc01=CCI=case, GSM4276357=dhsc03=CCI=case, GSM4276358=dhsc05=CCI=case
    # GSM4276359=dhsc13=SHAM=control, GSM4276360=dhsc15=SHAM=control
    rows = [
        ("GSM4276356", "case"),     # dhsc01 CCI
        ("GSM4276357", "case"),     # dhsc03 CCI
        ("GSM4276358", "case"),     # dhsc05 CCI
        ("GSM4276359", "control"),  # dhsc13 SHAM
        ("GSM4276360", "control"),  # dhsc15 SHAM
    ]
    pd.DataFrame(rows, columns=["sample_id", "group"]).to_csv(
        GEO_CACHE / "GSE143895_groups.csv", index=False)
    logger.info("  Restored GSM-ID groups.csv: 3 CCI case, 2 SHAM control")
    logger.info("  Title matching: 'dhsc01' ↔ 'dhsc01_abundance.tsv' will work")


# ---------------------------------------------------------------------------
# GSE113941 — CIPN vs Normal (input fraction, DRG)
# ---------------------------------------------------------------------------

def fix_gse113941() -> None:
    logger.info("\n=== GSE113941 (CIPN vs Normal, DRG input fraction) ===")
    src = GEO_CACHE / "GSE113941" / "GSE113941_uniform_process_tpm_vals.txt.gz"
    if not src.exists():
        logger.warning("  Source not found: %s", src)
        return

    df = pd.read_csv(src, sep="\t")
    # Keep gene_short_name + input fraction columns only
    # Normal input: columns containing "normal_input"
    # CIPN input: columns containing "CIPN_input" but NOT "subsample"
    keep_cols = ["gene_short_name"]
    rename_map = {}
    ctrl_count = 0
    case_count = 0
    for col in df.columns:
        if "normal_input" in col and "TRAP" not in col:
            ctrl_count += 1
            new_name = f"normal_control_{ctrl_count}"
            keep_cols.append(col)
            rename_map[col] = new_name
        elif ("CIPN_input" in col and "subsample" not in col
              and "MERGED" not in col and "TRAP" not in col):
            case_count += 1
            new_name = f"paclitaxel_case_{case_count}"
            keep_cols.append(col)
            rename_map[col] = new_name

    df_clean = df[keep_cols].rename(columns=rename_map)
    # Set gene_short_name as first column
    out = GEO_CACHE / "GSE113941" / "GSE113941_raw_count_input.txt.gz"
    save_gz(df_clean, out)
    # groups.csv uses GSM IDs (managed separately) — count file columns use vocab tokens
    logger.info("  Count file: %d control (normal), %d case (paclitaxel/CIPN)",
                ctrl_count, case_count)


# ---------------------------------------------------------------------------
# GSE306455 — SNI vs Sham spinal cord (exclude LDN treatment group)
# ---------------------------------------------------------------------------

def fix_gse306455() -> None:
    logger.info("\n=== GSE306455 (SNI vs Sham spinal cord) ===")
    src = GEO_CACHE / "GSE306455" / "GSE306455_allsamples.count.txt.gz"
    if not src.exists():
        logger.warning("  Source not found: %s", src)
        return

    df = pd.read_csv(src, sep="\t")
    # Columns: Geneid, L1, L2, L3 (LDN+SNI), S1, S2, S3 (Sham), N1, N2, N3 (SNI)
    # Keep only SNI (N) and Sham (S), rename to vocab-matchable names
    rename_map = {
        "S1": "sham_control_1", "S2": "sham_control_2", "S3": "sham_control_3",
        "N1": "sni_case_1",     "N2": "sni_case_2",     "N3": "sni_case_3",
    }
    cols_keep = ["Geneid", "S1", "S2", "S3", "N1", "N2", "N3"]
    df_clean = df[cols_keep].rename(columns=rename_map)
    out = GEO_CACHE / "GSE306455" / "GSE306455_raw_count_sni_sham.txt.gz"
    save_gz(df_clean, out)
    # groups.csv uses GSM IDs (managed separately) — count file columns use vocab tokens
    logger.info("  Count file: 3 sham_control, 3 sni_case (LDN excluded)")


# ---------------------------------------------------------------------------
# GSE289097 — CFA vs CON (trigeminal ganglion)
# ---------------------------------------------------------------------------

def fix_gse289097() -> None:
    logger.info("\n=== GSE289097 (CFA vs CON, trigeminal ganglion) ===")
    src = GEO_CACHE / "GSE289097" / "GSE289097_Expression_Gene.xlsx"
    if not src.exists():
        logger.warning("  Source not found: %s", src)
        return

    # Read with openpyxl — data starts at row 9 (0-indexed row 8)
    import openpyxl
    wb = openpyxl.load_workbook(src, read_only=True)
    ws = wb.active
    rows_all = list(ws.iter_rows(values_only=True))
    wb.close()

    header_idx = next(i for i, r in enumerate(rows_all) if r[0] == "Track_id")
    header = list(rows_all[header_idx])
    data   = [list(r) for r in rows_all[header_idx + 1:] if any(v is not None for v in r)]

    df = pd.DataFrame(data, columns=header)
    df = df.dropna(subset=["Gene_Name"])

    # Columns: Track_id, Gene_Name, Locus, Strand, Gene_Type, C1, C2, C3, N1, N2, N3
    # From GEO: C1/C2/C3 = CFA (case), N1/N2/N3 = CON/normal (control)
    # Rename to vocab-matchable names
    rename_map = {
        "C1": "cfa_case_1", "C2": "cfa_case_2", "C3": "cfa_case_3",
        "N1": "naive_control_1", "N2": "naive_control_2", "N3": "naive_control_3",
    }
    df_out = df[["Gene_Name", "C1", "C2", "C3", "N1", "N2", "N3"]].rename(columns=rename_map)
    out = GEO_CACHE / "GSE289097" / "GSE289097_raw_count_fpkm.txt.gz"
    save_gz(df_out, out)
    # groups.csv uses GSM IDs (managed separately) — count file columns use vocab tokens
    logger.info("  Count file: 3 cfa_case, 3 naive_control")


# ---------------------------------------------------------------------------
# GSE162284 — FLIT vs Sham (cortical neuropathic pain, Day3 vs Day0)
# ---------------------------------------------------------------------------

def fix_gse162284() -> None:
    logger.info("\n=== GSE162284 (FLIT vs Sham, cortical, Day3 vs Day0) ===")
    src = GEO_CACHE / "GSE162284" / "GSE162284_processed_RNAseq_rawcounts_12082022.csv.gz"
    if not src.exists():
        logger.warning("  Source not found: %s", src)
        return

    df = pd.read_csv(src)
    # Columns: Unnamed:0 (gene), Sham1_Day0..Sham4_Day0, FLIT1_Day3..FLIT4_Day3,
    #          FLIT1_Day21..FLIT4_Day21
    # Use Day3 FLIT (most acute neuropathic pain timepoint) vs Sham Day0
    gene_col = df.columns[0]
    sham_cols = [c for c in df.columns if "Sham" in c and "Day0" in c]
    flit_cols = [c for c in df.columns if "FLIT" in c and "Day3" in c]

    rename_map = {c: f"sham_control_{i+1}" for i, c in enumerate(sham_cols)}
    rename_map.update({c: f"injured_case_{i+1}" for i, c in enumerate(flit_cols)})

    df_out = df[[gene_col] + sham_cols + flit_cols].rename(columns=rename_map)
    df_out = df_out.rename(columns={gene_col: "gene_id"})
    out = GEO_CACHE / "GSE162284" / "GSE162284_raw_count_day3.txt.gz"
    save_gz(df_out, out)

    # groups.csv uses GSM IDs (managed separately) — count file columns use vocab tokens
    logger.info("  Count file: %d sham_control (Day0), %d injured_case (FLIT Day3)",
                len(sham_cols), len(flit_cols))


# ---------------------------------------------------------------------------
# GSE276193 — endometriosis vs non-endometriosis follicular fluid
# ---------------------------------------------------------------------------

def _series_samples(accession: str) -> list[tuple[str, str, str]]:
    """(GSM, title, supplementary filename) from the cached series matrix."""
    path = GEO_CACHE / f"{accession}_series_matrix.txt.gz"
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            cells = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            fields.setdefault(cells[0], [v.strip() for v in cells[1:]])
    return list(zip(fields["!Sample_geo_accession"], fields["!Sample_title"],
                    (Path(s).name for s in fields["!Sample_supplementary_file_1"])))


def fix_gse276193() -> None:
    """Rename the 18 sample columns to GSM ids and write groups.csv.

    The deposited matrix carries 14 annotation columns (gene_name, description,
    length, per-arm means, and seven functional-annotation columns) alongside
    18 sample columns named EM1..EM8 and CON1..CON10. Neither the sample titles
    ("endometriosis, follicle fluid, rep1") nor any vocabulary token matches
    those column names, so none of the DA's three group-assignment routes can
    reach the study as deposited.

    The mapping is not guessed: GEO's own per-sample supplementary filenames
    are `GSM8493022_EM1_processedfile.csv.gz`, so the column name is the middle
    token, and the arm is in the title. Deriving both from the series matrix
    rather than writing a table out by hand means the study cannot drift from
    what GEO says about it.
    """
    logger.info("\n=== GSE276193 (endometriosis vs non-endometriosis, follicular fluid) ===")
    src = GEO_CACHE / "GSE276193" / "GSE276193_merged_counts.csv.gz"
    if not src.exists():
        logger.warning("  Source not found: %s — run 05c first", src)
        return

    column_of, arm_of = {}, {}
    for gsm, title, filename in _series_samples("GSE276193"):
        parts = filename.split("_")
        if len(parts) < 2:
            logger.warning("  %s: cannot read a column name from %s", gsm, filename)
            return
        column_of[parts[1]] = gsm
        # "Non-endometriosis" contains "endometriosis", so test the negative first.
        arm_of[gsm] = "control" if title.lower().startswith("non-") else "case"

    frame = pd.read_csv(src, low_memory=False)
    gene_col = frame.columns[0]
    present = [c for c in frame.columns if c in column_of]
    missing = sorted(set(column_of) - set(present))
    if missing:
        # Silently pooling fewer samples than the arm counts claim is the
        # failure this whole file exists to prevent. Refuse instead.
        logger.warning("  columns absent from the matrix: %s", ", ".join(missing))
        return

    out = frame[[gene_col] + present].rename(columns={gene_col: "gene_id"})
    out = out.rename(columns={c: column_of[c] for c in present})
    save_gz(out, GEO_CACHE / "GSE276193" / "GSE276193_raw_count_merged.txt.gz")

    rows = [(gsm, arm_of[gsm]) for gsm in out.columns[1:]]
    pd.DataFrame(rows, columns=["sample_id", "group"]).to_csv(
        GEO_CACHE / "GSE276193_groups.csv", index=False)
    n_case = sum(1 for _g, a in rows if a == "case")
    logger.info("  groups.csv: %d case (endometriosis), %d control",
                n_case, len(rows) - n_case)

    # The source has to go. Both filenames match find_suppl_count_file()'s raw
    # tier, and that tier deliberately has no preference by filename -- for a
    # good reason, recorded there -- so leaving both in place makes the choice
    # alphabetical, which picks the deposit and then fails to assign groups
    # because its columns are EM1/CON1. Step 05c will fetch it again on demand
    # and this step will consume it again; what must not persist is two files in
    # one tier where only one can be analysed.
    src.unlink()
    logger.info("  Removed the source deposit: both names match the raw tier, "
                "where file choice is alphabetical")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    fix_gse143895()
    fix_gse113941()
    fix_gse306455()
    fix_gse289097()
    fix_gse162284()
    fix_gse276193()
    logger.info("\nPreprocessing complete.")


if __name__ == "__main__":
    main()
