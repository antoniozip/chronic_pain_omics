"""Accession-keyed resolver for pending screening records (Task 2 of
plan/2026-10-30-resolve-pending-screening.md)."""
import subprocess
from pathlib import Path


def _run(decisions: Path, sheet: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["python", "scripts/resolve_pending_candidates.py",
         "--decisions", str(decisions), "--sheet", str(sheet)],
        capture_output=True, text=True, timeout=120)


def test_unknown_accession_rejected(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text("accession,verdict,reason_code\nGSE1,pending,usable_pending_file_check\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE999,exclude,no_count_matrix_in_deposit\n")
    res = _run(dec, sheet)
    assert res.returncode != 0
    assert "unknown" in (res.stdout + res.stderr)


def test_non_pending_row_refused(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text("accession,verdict,reason_code\nGSE1,include,\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE1,exclude,no_count_matrix_in_deposit\n")
    res = _run(dec, sheet)
    assert res.returncode != 0
    assert sheet.read_text() == "accession,verdict,reason_code\nGSE1,include,\n"


def test_exclude_requires_reason_code(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text("accession,verdict,reason_code\nGSE1,pending,usable_pending_file_check\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE1,exclude,\n")
    res = _run(dec, sheet)
    assert res.returncode != 0


def test_include_pending_row_updated(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text(
        "accession,verdict,reason_code\nGSE1,pending,usable_pending_file_check\n"
        "GSE2,exclude,no_count_matrix_in_deposit\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE1,include,\n")
    res = _run(dec, sheet)
    assert res.returncode == 0
    assert sheet.read_text() == ("accession,verdict,reason_code\n"
                                 "GSE1,include,\nGSE2,exclude,no_count_matrix_in_deposit\n")


def test_include_to_exclude_requires_flag(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text("accession,verdict,reason_code\nGSE1,include,\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE1,exclude,step05_error_row_id_length_mismatch\n")
    res = _run(dec, sheet)
    assert res.returncode != 0
    assert sheet.read_text() == "accession,verdict,reason_code\nGSE1,include,\n"


def test_include_to_exclude_with_flag(tmp_path):
    sheet = tmp_path / "candidates.csv"
    sheet.write_text("accession,verdict,reason_code\nGSE1,include,\n")
    dec = tmp_path / "d.csv"
    dec.write_text("accession,verdict,reason_code\nGSE1,exclude,step05_error_row_id_length_mismatch\n")
    res = subprocess.run(
        ["python", "scripts/resolve_pending_candidates.py",
         "--decisions", str(dec), "--sheet", str(sheet), "--allow-include-to-exclude"],
        capture_output=True, text=True, timeout=120)
    assert res.returncode == 0
    assert sheet.read_text() == ("accession,verdict,reason_code\n"
                                 "GSE1,exclude,step05_error_row_id_length_mismatch\n")
