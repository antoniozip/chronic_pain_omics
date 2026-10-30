"""The no-split review interface: what the server answers, and what it refuses to save.

The interface exists because `no_case_control_split_in_series_matrix` cannot tell
a deposit with no contrast from a parser that could not read one. These tests
use a miniature repository so they need neither the GEO cache nor the sheet.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from tests.conftest import load_script_module

rv = load_script_module("review_no_split.py")

MATRIX = (
    '!Series_title\t"A study"\n'
    '!Series_summary\t"Dorsal root ganglia were profiled."\n'
    '!Series_overall_design\t"Four treated, four DMSO."\n'
    '!Sample_title\t"s1"\t"s2"\t"s3"\t"s4"\t"s5"\t"s6"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\t"GSM5"\t"GSM6"\n'
    '!Sample_source_name_ch1\t"DRG"\t"DRG"\t"DRG"\t"DRG"\t"DRG"\t"DRG"\n'
    '!Sample_characteristics_ch1\t"treatment: drugX"\t"treatment: drugX"'
    '\t"treatment: drugX"\t"treatment: DMSO"\t"treatment: DMSO"\t"treatment: DMSO"\n'
    '!Sample_characteristics_ch1\t"site: L4"\t"site: L5"\t"site: L4"\t"site: L5"'
    '\t"site: L4"\t"site: L5"\n'
    "!series_matrix_table_begin\n"
)


@pytest.fixture()
def ctx(tmp_path: Path):
    prisma = tmp_path / "literature" / "prisma"
    prisma.mkdir(parents=True)
    cols = ["accession", "species", "n_samples", "year", "title", "url",
            "verdict", "reason_code"]
    rows = [
        ["GSE1", "Mus musculus", 6, 2020, "A study", "https://x/1", "exclude", rv.REASON],
        ["GSE2", "Homo sapiens", 9, 2021, "Other", "https://x/2", "include", "fine"],
        ["GSE3", "Mus musculus", 6, 2022, "Third", "https://x/3", "exclude", rv.REASON],
    ]
    with open(tmp_path / rv.SHEET, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows)
    with open(tmp_path / rv.AUDIT, "w", newline="") as fh:
        fh.write("accession,status,field,n_case,n_control,case_label,control_label,"
                 "field_looks_non_disease\nGSE1,candidate_contrast_found,treatment,3,3,"
                 "drugX,DMSO,no\n")
    cache = tmp_path / rv.CACHE
    cache.mkdir(parents=True)
    for acc in ("GSE1", "GSE3"):
        with gzip.open(cache / f"{acc}_series_matrix.txt.gz", "wt") as fh:
            fh.write(MATRIX)
    rv._load.cache_clear()
    reviews = tmp_path / "reviews.csv"
    return rv.Context(tmp_path, reviews, rv.build_queue(tmp_path, reviews))


RULE = {"case": {"values": ["drugX"]}, "control": {"values": ["DMSO"]}}


def _field(ctx, label: str) -> str:
    return next(f["id"] for f in rv.series_detail(ctx, "GSE1")["fields"]
                if f["label"] == label)


def test_queue_holds_only_the_series_excluded_for_this_reason(ctx):
    assert [q["accession"] for q in ctx.queue] == ["GSE1", "GSE3"]


def test_queue_carries_the_audit_and_marks_a_series_never_audited(ctx):
    by = {q["accession"]: q for q in ctx.queue}
    assert by["GSE1"]["audit_status"] == "candidate_contrast_found"
    assert by["GSE3"]["audit_status"] == "not_audited"


def test_detail_shows_the_deposits_own_text_and_every_field(ctx):
    d = rv.series_detail(ctx, "GSE1")
    assert d["header"]["Series_overall_design"] == "Four treated, four DMSO."
    assert [f["label"] for f in d["fields"]] == [
        "title", "source_name", "char: treatment", "char: site"]
    assert d["gsm"][0] == "GSM1"


def test_a_field_the_parser_skips_says_which_word_skipped_it(ctx):
    site = next(f for f in rv.series_detail(ctx, "GSE1")["fields"]
                if f["label"] == "char: site")
    assert site["skipped_by"] == ["site"]


def test_a_field_without_variation_is_flagged_constant(ctx):
    src = next(f for f in rv.series_detail(ctx, "GSE1")["fields"]
               if f["label"] == "source_name")
    assert src["constant"]


def test_rule_counts_exact_values_and_regexes_alike(ctx):
    values = ["drugX", "drugX", "drugX", "DMSO", "DMSO", "DMSO"]
    exact = rv.apply_rule(values, RULE)
    assert (exact["n_case"], exact["n_control"]) == (3, 3)
    rx = rv.apply_rule(values, {"case": {"regex": "^DRUG"}, "control": {"regex": "dms"}})
    assert (rx["n_case"], rx["n_control"]) == (3, 3)


def test_a_sample_matching_both_arms_is_a_conflict_in_neither(ctx):
    res = rv.apply_rule(["a b"], {"case": {"regex": "a"}, "control": {"regex": "b"}})
    assert (res["n_case"], res["n_control"], res["n_conflict"]) == (0, 0, 1)


def test_an_invalid_regex_is_refused_not_ignored(ctx):
    with pytest.raises(ValueError, match="invalid regex"):
        rv.preview(ctx, "GSE1", _field(ctx, "char: treatment"), {"case": {"regex": "("}})


def test_saving_a_contrast_records_counts_the_server_computed(ctx):
    row = rv.save_review(ctx, {"accession": "GSE1", "decision": "contrast_exists",
                               "field": _field(ctx, "char: treatment"), "rule": RULE,
                               "note": "treated vs vehicle"})
    assert (row["n_case"], row["n_control"], row["field"]) == ("3", "3", "char: treatment")
    saved = list(csv.DictReader(open(ctx.reviews)))
    assert saved[0]["accession"] == "GSE1"
    assert json.loads(saved[0]["rule_json"])["case"]["values"] == ["drugX"]


def test_a_contrast_with_one_sample_in_an_arm_is_refused(ctx):
    rule = {"case": {"values": ["drugX"]}, "control": {"regex": "^DMSO$"}}
    thin = {"case": {"regex": "^drugX$"}, "control": {"values": ["nonexistent"]}}
    for r in (thin, {"case": rule["case"], "control": {}}):
        with pytest.raises(ValueError, match="at least 2"):
            rv.save_review(ctx, {"accession": "GSE1", "decision": "contrast_exists",
                                 "field": _field(ctx, "char: treatment"), "rule": r})
    assert not ctx.reviews.exists()


def test_a_contrast_that_overlaps_itself_is_refused(ctx):
    both = {"case": {"regex": "."}, "control": {"regex": "."}}
    with pytest.raises(ValueError, match="both case and control"):
        rv.save_review(ctx, {"accession": "GSE1", "decision": "contrast_exists",
                             "field": _field(ctx, "char: treatment"), "rule": both})


def test_no_contrast_saves_without_a_rule(ctx):
    row = rv.save_review(ctx, {"accession": "GSE3", "decision": "no_contrast",
                               "note": "genotype only"})
    assert row["rule_json"] == "" and row["n_case"] == ""


def test_saving_twice_replaces_the_decision(ctx):
    rv.save_review(ctx, {"accession": "GSE3", "decision": "unclear"})
    rv.save_review(ctx, {"accession": "GSE3", "decision": "no_contrast"})
    saved = list(csv.DictReader(open(ctx.reviews)))
    assert [(r["accession"], r["decision"]) for r in saved] == [("GSE3", "no_contrast")]


def test_an_accession_outside_the_queue_is_refused(ctx):
    with pytest.raises(KeyError):
        rv.save_review(ctx, {"accession": "GSE2", "decision": "no_contrast"})
    with pytest.raises(KeyError):
        rv.series_detail(ctx, "../../etc/passwd")


def test_an_unknown_decision_is_refused(ctx):
    with pytest.raises(ValueError, match="decision must be"):
        rv.save_review(ctx, {"accession": "GSE1", "decision": "maybe"})


def test_the_review_never_touches_the_screening_sheet(ctx):
    sheet = ctx.root / rv.SHEET
    before = sheet.read_bytes()
    rv.save_review(ctx, {"accession": "GSE1", "decision": "no_contrast"})
    assert sheet.read_bytes() == before


def test_a_saved_review_comes_back_in_the_queue(ctx):
    rv.save_review(ctx, {"accession": "GSE1", "decision": "unclear", "note": "n"})
    by = {q["accession"]: q for q in rv.build_queue(ctx.root, ctx.reviews)}
    assert by["GSE1"]["review"]["decision"] == "unclear"
    assert by["GSE3"]["review"] is None
