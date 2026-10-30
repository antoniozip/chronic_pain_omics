"""Applying the count-column review: what reaches the cache and the sheet."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from tests.conftest import load_script_module

mc = load_script_module("map_count_columns.py")
ap = load_script_module("apply_count_mapping.py")
rv = mc.rv

MATRIX = (
    '!Series_title\t"A study"\n'
    '!Sample_title\t"c1"\t"c2"\t"k1"\t"k2"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\n'
    "!series_matrix_table_begin\n"
)
COLS = ["accession", "species", "n_samples", "year", "title", "url", "verdict", "reason_code",
        "reason_source", "annotation_status", "n_case", "n_control"]


@pytest.fixture()
def root(tmp_path: Path) -> Path:
    (tmp_path / "literature" / "prisma").mkdir(parents=True)
    rows = [["GSE1", "Mus musculus", 4, 2021, "A", "u", "pending", mc.REASON, "x", "", "", ""],
            ["GSE2", "Mus musculus", 4, 2021, "B", "u", "pending", mc.REASON, "x", "", "", ""],
            ["GSE3", "Mus musculus", 4, 2021, "C", "u", "include", "fine", "x", "", "", ""]]
    with open(tmp_path / rv.SHEET, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\r\n")
        w.writerow(COLS)
        w.writerows(rows)
    cache = tmp_path / rv.CACHE
    for acc in ("GSE1", "GSE2"):
        (cache / acc).mkdir(parents=True)
        with gzip.open(cache / f"{acc}_series_matrix.txt.gz", "wt") as fh:
            fh.write(MATRIX)
        (cache / f"{acc}_groups.csv").write_text(
            "sample_id,group\nGSM1,control\nGSM2,control\nGSM3,case\nGSM4,case\n")
    pd.DataFrame({"sym": ["A", "B", "B"], "k1": [1.0, 4.0, 9.0], "k2": [2.0, 4.0, 9.0],
                  "c1": [3.0, 5.0, 1.0], "c2": [4.0, 6.0, 1.0], "extra": [0, 0, 0]}).to_csv(
        cache / "GSE1" / "GSE1_counts.txt.gz", sep="\t", index=False)
    pd.DataFrame({"sym": ["A"], "k1": [1], "k2": [1], "c1": [1], "c2": [1]}).to_csv(
        cache / "GSE2" / "GSE2_counts.txt.gz", sep="\t", index=False)
    rv._load.cache_clear()
    return tmp_path


def review(root: Path, *rows: dict[str, str]) -> Path:
    path = root / "reviews.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=mc.COLUMNS)
        w.writeheader()
        w.writerows({k: "" for k in mc.COLUMNS} | r for r in rows)
    return path


OK = {"accession": "GSE1", "decision": "mapping_ok", "gene_col": "sym", "scale": "counts",
      "mapping_json": json.dumps({"c1": "GSM1", "c2": "GSM2", "k1": "GSM3", "k2": "GSM4"}),
      "drop_json": "[]", "n_case": "2", "n_control": "2"}


def sheet_row(root: Path, acc: str) -> dict[str, str]:
    raw = (root / rv.SHEET).read_bytes().decode()
    return next(r for r in csv.DictReader(raw.splitlines()) if r["accession"] == acc)


def test_a_dry_run_changes_nothing(root: Path):
    before = (root / rv.SHEET).read_bytes()
    done = ap.apply(root, review(root, OK), write=False)
    assert done["mapped"] == ["GSE1"]
    assert (root / rv.SHEET).read_bytes() == before
    assert (root / rv.CACHE / "GSE1" / "GSE1_counts.txt.gz").exists()


def test_the_table_is_renamed_to_gsm_ids_and_the_deposit_moved_aside(root: Path):
    ap.apply(root, review(root, OK), write=True)
    out = pd.read_csv(root / rv.CACHE / "GSE1" / "GSE1_raw_count_merged.txt.gz", sep="\t")
    assert list(out.columns) == ["gene_id", "GSM1", "GSM2", "GSM3", "GSM4"]
    assert not (root / rv.CACHE / "GSE1" / "GSE1_counts.txt.gz").exists()
    assert (root / ap.ORIGINALS / "GSE1" / "GSE1_counts.txt.gz").exists()


def test_a_duplicated_gene_label_keeps_its_highest_mean_row(root: Path):
    ap.apply(root, review(root, OK), write=True)
    out = pd.read_csv(root / rv.CACHE / "GSE1" / "GSE1_raw_count_merged.txt.gz", sep="\t")
    assert out["gene_id"].tolist() == ["A", "B"]
    assert out.loc[out["gene_id"] == "B", "GSM3"].item() == 9.0


def test_log2_values_come_back_to_the_linear_scale():
    df = pd.DataFrame({"g": ["a"], "x": [3.0], "y": [0.0]})
    rev = {"mapping_json": json.dumps({"x": "GSM1", "y": "GSM2"}), "drop_json": "[]",
           "scale": "log2", "gene_col": "g"}
    out = ap.build_table(df, rev)
    assert out["GSM1"].item() == 7.0 and out["GSM2"].item() == 0.0


def test_the_sheet_row_goes_to_pending_with_the_reviewers_counts_and_keeps_crlf(root: Path):
    ap.apply(root, review(root, OK), write=True)
    row = sheet_row(root, "GSE1")
    assert (row["verdict"], row["reason_code"]) == ("pending", ap.PENDING)
    assert (row["n_case"], row["n_control"]) == ("2", "2")
    assert b"\r\n" in (root / rv.SHEET).read_bytes()
    assert sheet_row(root, "GSE3")["verdict"] == "include"


def test_dropped_samples_leave_the_groups_file_and_the_original_is_kept(root: Path):
    rev = OK | {"mapping_json": json.dumps({"c1": "GSM1", "c2": "GSM2", "k1": "GSM3"}),
                "drop_json": json.dumps(["GSM4"]), "n_case": "1", "n_control": "2"}
    ap.apply(root, review(root, rev), write=True)
    groups = (root / rv.CACHE / "GSE1_groups.csv").read_text()
    assert "GSM4" not in groups and "GSM3" in groups
    assert "GSM4" in (root / ap.ORIGINALS / "GSE1" / "GSE1_groups.csv").read_text()


def test_an_ineligible_series_is_excluded_without_touching_its_data(root: Path):
    rev = {"accession": "GSE2", "decision": "ineligible",
           "reason_code": "in_vitro_stimulus_not_patient_contrast", "note": "organ culture"}
    done = ap.apply(root, review(root, rev), write=True)
    assert done["excluded"] == ["GSE2"]
    row = sheet_row(root, "GSE2")
    assert (row["verdict"], row["reason_code"]) == (
        "exclude", "in_vitro_stimulus_not_patient_contrast")
    assert (root / rv.CACHE / "GSE2" / "GSE2_counts.txt.gz").exists()


def test_unclear_changes_nothing(root: Path):
    before = (root / rv.SHEET).read_bytes()
    done = ap.apply(root, review(root, {"accession": "GSE1", "decision": "unclear"}), True)
    assert done == {"mapped": [], "excluded": []}
    assert (root / rv.SHEET).read_bytes() == before


def test_a_series_no_longer_pending_for_this_reason_is_refused(root: Path):
    with pytest.raises(ValueError, match="refusing to overwrite"):
        ap.apply(root, review(root, OK | {"accession": "GSE3"}), write=True)


def test_counts_that_disagree_with_the_review_are_refused(root: Path):
    with pytest.raises(ValueError, match="columns built"):
        ap.apply(root, review(root, OK | {"n_case": "3"}), write=True)
