"""The count-column mapping interface: what it proposes, and what it refuses to save.

A wrong pairing produces an effect table of exactly the right shape with the
wrong samples in each arm, so these tests pin the cases where a proposal must
stay silent as firmly as the cases where it must speak. They use a miniature
repository, so they need neither the GEO cache nor the screening sheet.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from tests.conftest import load_script_module

mc = load_script_module("map_count_columns.py")
rv = mc.rv

MATRIX = (
    '!Series_title\t"A study"\n'
    '!Series_overall_design\t"Four animals per arm."\n'
    '!Sample_title\t"DRG-Con_rep1"\t"DRG-Con_rep2"\t"DRG-CPIP_rep1"\t"DRG-CPIP_rep2"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\n'
    '!Sample_source_name_ch1\t"DRG"\t"DRG"\t"DRG"\t"DRG"\n'
    "!series_matrix_table_begin\n"
)


def samples(*titles: str) -> list[dict]:
    return [{"gsm": f"GSM{i + 1}", "strings": [t]} for i, t in enumerate(titles)]


@pytest.fixture()
def ctx(tmp_path: Path):
    prisma = tmp_path / "literature" / "prisma"
    prisma.mkdir(parents=True)
    cols = ["accession", "species", "n_samples", "year", "title", "url", "verdict",
            "reason_code"]
    with open(tmp_path / rv.SHEET, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerow(["GSE1", "Rattus norvegicus", 4, 2021, "A study", "https://x/1",
                    "pending", mc.REASON])
        w.writerow(["GSE2", "Mus musculus", 4, 2022, "Other", "https://x/2", "include", "ok"])
    cache = tmp_path / rv.CACHE
    (cache / "GSE1").mkdir(parents=True)
    with gzip.open(cache / "GSE1_series_matrix.txt.gz", "wt") as fh:
        fh.write(MATRIX)
    (cache / "GSE1_groups.csv").write_text(
        "sample_id,group\nGSM1,control\nGSM2,control\nGSM3,case\nGSM4,case\n")
    pd.DataFrame({"Gene Symbol": ["A", "B"], "Con1 FPKM": [1.5, 2.5], "Con2 FPKM": [1.0, 3.0],
                  "CPIP1 FPKM": [5.5, 1.0], "CPIP2 FPKM": [6.5, 0.5]}).to_csv(
        cache / "GSE1" / "GSE1_Raw_gene_FPKM.txt.gz", sep="\t", index=False)
    rv._load.cache_clear()
    reviews = tmp_path / "reviews.csv"
    return mc.Context(tmp_path, reviews, mc.build_queue(tmp_path, reviews))


GOOD = {"Con1 FPKM": "GSM1", "Con2 FPKM": "GSM2", "CPIP1 FPKM": "GSM3", "CPIP2 FPKM": "GSM4"}


def rec(**kw):
    base = {"accession": "GSE1", "decision": "mapping_ok", "gene_col": "Gene Symbol",
            "scale": "linear", "mapping": GOOD}
    return {**base, **kw}


def test_label_key_drops_the_words_that_decorate_a_label():
    assert mc.label_key("CPIP1 FPKM") == "cpip1"
    assert mc.label_key("DRG-CPIP_rep1") == "drgcpip1"
    assert mc.label_key("norm.CT1") == "ct1"


def test_equal_labels_pair_whatever_the_punctuation():
    got = mc.propose(["NE_1", "NE_2"], samples("NE-2", "NE-1"))
    assert got["NE_1"] == {"gsm": "GSM2", "how": "equal"}
    assert got["NE_2"] == {"gsm": "GSM1", "how": "equal"}


def test_a_label_contained_in_a_title_pairs():
    got = mc.propose(["CPIP1 FPKM", "Con1 FPKM"], samples("DRG-Con_rep1", "DRG-CPIP_rep1"))
    assert (got["CPIP1 FPKM"]["gsm"], got["Con1 FPKM"]["gsm"]) == ("GSM2", "GSM1")
    assert got["Con1 FPKM"]["how"] == "contained"


def test_a_subsequence_pairs_only_when_the_numbers_agree():
    got = mc.propose(["norm.CT1", "norm.EN2"], samples("CTR1", "END2"))
    assert got["norm.CT1"]["how"] == "subsequence"
    other = mc.propose(["norm.CT1"], samples("CTR2"))
    assert other["norm.CT1"]["how"] == "positional"   # by name it would not have paired


def test_a_column_that_two_samples_claim_equally_is_not_paired_by_name():
    # "c1" sits inside both titles: the name cannot decide, so only position can
    got = mc.propose(["C1"], samples("C1_PVT", "C1_pIC"))
    assert got == {}


def test_leftovers_are_paired_by_position_and_marked_as_such():
    got = mc.propose(["a", "b"], samples("zz", "yy"))
    assert [v["how"] for v in got.values()] == ["positional", "positional"]
    assert [v["gsm"] for v in got.values()] == ["GSM1", "GSM2"]


def test_positional_pairing_takes_the_raw_count_block_among_several():
    cols = ([f"fpkm_{i}" for i in range(2)] + [f"read_count_{i}" for i in range(2)]
            + [f"norm_read_count_{i}" for i in range(2)])
    got = mc.propose(cols, samples("zz", "yy"))
    assert sorted(got) == ["read_count_0", "read_count_1"]


def test_positional_pairing_refuses_blocks_it_cannot_tell_apart():
    cols = [f"fpkm_{i}" for i in range(2)] + [f"tpm_{i}" for i in range(2)]
    assert mc.propose(cols, samples("zz", "yy")) == {}


def test_the_05d_log_gives_the_gsm_behind_each_group_named_column(tmp_path: Path):
    log = tmp_path / "run_05d.log"
    log.write_text("=== GSE9 ===\n  GSM11 → oa_case_1 (100 genes)\n"
                   "  GSM12 → vehicle_control_1 (100 genes)\n=== GSE10 ===\n"
                   "  GSM99 → oa_case_1 (5 genes)\n")
    assert mc.read_05d_log(log, "GSE9") == {"oa_case_1": "GSM11", "vehicle_control_1": "GSM12"}
    assert mc.read_05d_log(log, "GSE404") == {}


def test_a_table_headed_on_its_second_row_is_reheaded(tmp_path: Path):
    path = tmp_path / "t.txt.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("Gene symbol\tNormalized RC (log2)\tUnnamed: 2\n\tWT_1\tWT_2\nXkr4\t6.1\t6.5\n")
    df = mc.read_table(path)
    assert list(df.columns) == ["Gene symbol", "WT_1", "WT_2"]
    assert len(df) == 1


def test_scale_is_read_off_the_columns_that_will_be_used():
    df = pd.DataFrame({"g": ["a", "b"], "x": [3, 9], "start": [10**8, 10**9]})
    assert mc.guess_scale(df, ["x"])["scale_guess"] == "counts"
    assert mc.guess_scale(df.assign(x=[1.5, 9.2]), ["x"])["scale_guess"] == "log2"
    assert mc.guess_scale(df.assign(x=[1.5, 900.2]), ["x"])["scale_guess"] == "linear"


def test_queue_holds_only_the_series_pending_for_this_reason(ctx):
    assert [q["accession"] for q in ctx.queue] == ["GSE1"]


def test_detail_proposes_the_pairing_and_shows_the_group(ctx):
    d = mc.study_detail(ctx, "GSE1")
    assert d["proposal"]["CPIP2 FPKM"]["gsm"] == "GSM4"
    assert {s["gsm"]: s["group"] for s in d["samples"]}["GSM1"] == "control"


def test_saving_a_complete_mapping_records_the_counts_the_server_computed(ctx):
    row = mc.save_review(ctx, rec())
    assert (row["n_case"], row["n_control"]) == ("2", "2")
    assert json.loads(next(csv.DictReader(open(ctx.reviews)))["mapping_json"]) == GOOD


def test_two_columns_for_one_sample_are_refused(ctx):
    with pytest.raises(ValueError, match="same sample"):
        mc.save_review(ctx, rec(mapping={**GOOD, "Con2 FPKM": "GSM1"}))


def test_a_case_control_sample_without_a_column_is_refused_unless_dropped(ctx):
    short = {k: v for k, v in GOOD.items() if v != "GSM4"}
    with pytest.raises(ValueError, match="no column"):
        mc.save_review(ctx, rec(mapping=short))
    with pytest.raises(ValueError, match="at least 2"):
        mc.save_review(ctx, rec(mapping=short, drop=["GSM4"]))


def test_the_gene_column_cannot_be_a_sample_column(ctx):
    with pytest.raises(ValueError, match="gene column"):
        mc.save_review(ctx, rec(gene_col="Con1 FPKM"))


def test_an_unknown_sample_or_column_is_refused(ctx):
    with pytest.raises(ValueError, match="not samples"):
        mc.save_review(ctx, rec(mapping={**GOOD, "Con1 FPKM": "GSM404"}))
    with pytest.raises(ValueError, match="numeric columns"):
        mc.save_review(ctx, rec(mapping={**GOOD, "Nope": "GSM1"}))


def test_an_ineligible_series_needs_a_reason_and_a_note(ctx):
    with pytest.raises(ValueError, match="reason code"):
        mc.save_review(ctx, {"accession": "GSE1", "decision": "ineligible", "note": "x"})
    with pytest.raises(ValueError, match="note"):
        mc.save_review(ctx, {"accession": "GSE1", "decision": "ineligible",
                             "reason_code": "no_pain_phenotype"})
    row = mc.save_review(ctx, {"accession": "GSE1", "decision": "ineligible",
                               "reason_code": "no_pain_phenotype", "note": "diet only"})
    assert row["n_case"] == ""


def test_an_accession_outside_the_queue_is_refused(ctx):
    with pytest.raises(KeyError):
        mc.save_review(ctx, rec(accession="GSE2"))


def test_the_tool_never_touches_the_sheet_or_the_cache(ctx):
    sheet = (ctx.root / rv.SHEET).read_bytes()
    groups = (ctx.root / rv.CACHE / "GSE1_groups.csv").read_bytes()
    mc.save_review(ctx, rec())
    assert (ctx.root / rv.SHEET).read_bytes() == sheet
    assert (ctx.root / rv.CACHE / "GSE1_groups.csv").read_bytes() == groups
