"""Tests for pipeline/02_screen.py (no filesystem or network calls)."""

import csv as csv_mod

from cp_multiomics.vocabulary import PainVocabulary, load_pain_vocabulary


def _make_hit(uid: str, doi: str = "", source: str = "pubmed") -> dict:
    return {
        "uid": uid, "source": source, "doi": doi, "title": f"Title {uid}",
        "authors": "Author A", "journal": "J Pain", "year": "2023",
        "abstract": "Abstract text.", "query_label": "test",
        "retrieved_at": "2024-01-01T00:00:00+00:00",
    }


def test_dedup_removes_cross_run_doi_duplicates(pipeline_02):
    m = pipeline_02
    r1 = m.ScreeningRecord.from_hit(_make_hit("111", doi="10.1/abc"))
    r2 = m.ScreeningRecord.from_hit(_make_hit("111", doi="10.1/abc"))
    unique, n = m.deduplicate([r1, r2])
    assert n == 1 and len(unique) == 1


def test_dedup_uid_fallback(pipeline_02):
    m = pipeline_02
    r1 = m.ScreeningRecord.from_hit(_make_hit("111"))
    r2 = m.ScreeningRecord.from_hit(_make_hit("111"))
    unique, n = m.deduplicate([r1, r2])
    assert n == 1


def test_from_hit_preserves_all_fields(pipeline_02):
    m = pipeline_02
    hit = _make_hit("999", doi="10.1/xyz", source="europe_pmc")
    r = m.ScreeningRecord.from_hit(hit)
    assert r.uid == "999"
    assert r.doi == "10.1/xyz"
    assert r.source == "europe_pmc"
    assert r.screening_decision == ""


def test_write_and_load_roundtrip(pipeline_02, tmp_path):
    m = pipeline_02
    records = [m.ScreeningRecord.from_hit(_make_hit(str(i), doi=f"10.1/{i}")) for i in range(3)]
    records[0].screening_decision = "include"
    records[1].screening_decision = "exclude"
    records[1].exclusion_reason = "wrong_condition"

    csv_path = tmp_path / "screening.csv"
    m.write_screening_csv(records, csv_path)
    loaded = m.load_existing_screening(csv_path)

    assert len(loaded) == 3
    assert loaded["doi:10.1/0"].screening_decision == "include"
    assert loaded["doi:10.1/1"].exclusion_reason == "wrong_condition"
    assert loaded["doi:10.1/2"].screening_decision == ""


def test_apply_decisions_merges_correctly(pipeline_02, tmp_path):
    m = pipeline_02
    records = [m.ScreeningRecord.from_hit(_make_hit(str(i), doi=f"10.1/{i}")) for i in range(2)]
    screening_csv = tmp_path / "screening.csv"
    m.write_screening_csv(records, screening_csv)

    decisions_csv = tmp_path / "decisions.csv"
    with open(decisions_csv, "w", newline="") as f:
        fieldnames = ["uid", "source", "doi", "screening_decision",
                      "exclusion_reason", "screened_by", "screening_notes"]
        w = csv_mod.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerow({"uid": "0", "source": "pubmed", "doi": "10.1/0",
                    "screening_decision": "include", "exclusion_reason": "",
                    "screened_by": "AGZ", "screening_notes": ""})
        w.writerow({"uid": "1", "source": "pubmed", "doi": "10.1/1",
                    "screening_decision": "exclude", "exclusion_reason": "no_omics_data",
                    "screened_by": "AGZ", "screening_notes": "no raw data"})

    m.apply_decisions(screening_csv, decisions_csv)
    result = m.load_existing_screening(screening_csv)

    assert result["doi:10.1/0"].screening_decision == "include"
    assert result["doi:10.1/1"].screening_decision == "exclude"
    assert result["doi:10.1/1"].exclusion_reason == "no_omics_data"
    assert result["doi:10.1/1"].screened_by == "AGZ"


# --- appended: auto-screening (Task 8) -------------------------------------
VOCAB = PainVocabulary(core=("chronic pain",), conditions=("fibromyalgia",),
                       mechanisms=("allodynia",))
OMICS = ("rna-seq", "gwas", "proteomics")


def _rec(pipeline_02, title="", abstract=""):
    return pipeline_02.ScreeningRecord(
        uid="1", source="pubmed", doi="", title=title, authors="", journal="",
        year="2020", abstract=abstract, query_label="q", retrieved_at="",
    )


def test_record_with_pain_and_omics_is_left_for_a_human(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "RNA-seq of chronic pain in rats"), VOCAB, OMICS)
    assert decision == ""
    assert reason == ""


def test_record_without_any_omics_term_is_auto_excluded(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "Acupuncture for chronic pain: a trial"), VOCAB, OMICS)
    assert decision == "exclude"
    assert reason == "no_omics_term"


def test_record_without_any_pain_term_is_auto_excluded(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "RNA-seq of maize leaf under salt stress"), VOCAB, OMICS)
    assert decision == "exclude"
    assert reason == "not_pain_phenotype"


def test_auto_screen_never_auto_includes(pipeline_02):
    """PRISMA: inclusion is always a human decision."""
    for title in ("RNA-seq of chronic pain", "GWAS of fibromyalgia and allodynia"):
        decision, _ = pipeline_02.auto_screen(_rec(pipeline_02, title), VOCAB, OMICS)
        assert decision != "include"


def test_abstract_is_searched_not_just_the_title(pipeline_02):
    decision, _ = pipeline_02.auto_screen(
        _rec(pipeline_02, "A cohort study", "We performed proteomics in fibromyalgia patients"),
        VOCAB, OMICS)
    assert decision == ""


def test_matching_is_case_insensitive(pipeline_02):
    decision, _ = pipeline_02.auto_screen(_rec(pipeline_02, "GWAS OF CHRONIC PAIN"), VOCAB, OMICS)
    assert decision == ""


def test_existing_human_decisions_are_never_overwritten(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    kept = _rec(pipeline_02, "Acupuncture for chronic pain")   # would be auto-excluded
    kept.screening_decision = "include"
    kept.screened_by = "AG Zippo"
    pipeline_02.write_screening_csv([kept], csv_path)

    counts = pipeline_02.apply_auto_exclusions(csv_path, VOCAB, OMICS)

    records = list(pipeline_02.load_existing_screening(csv_path).values())
    assert records[0].screening_decision == "include"
    assert records[0].screened_by == "AG Zippo"
    assert counts["skipped_already_decided"] == 1


def test_auto_excluded_rows_are_attributed_and_reasoned(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    pipeline_02.write_screening_csv([_rec(pipeline_02, "maize leaf proteomics")], csv_path)

    pipeline_02.apply_auto_exclusions(csv_path, VOCAB, OMICS)

    rec = list(pipeline_02.load_existing_screening(csv_path).values())[0]
    assert rec.screening_decision == "exclude"
    assert rec.exclusion_reason == "not_pain_phenotype"
    assert rec.screened_by == "auto"


def test_new_reason_codes_are_registered(pipeline_02):
    assert "no_omics_term" in pipeline_02.EXCLUSION_REASONS
    assert "not_pain_phenotype" in pipeline_02.EXCLUSION_REASONS


def test_audit_sample_is_deterministic_and_only_auto_rows(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    recs = []
    for i in range(20):
        r = _rec(pipeline_02, f"maize study {i}")
        r.uid = str(i)          # dedup_key is source:uid when doi is empty
        recs.append(r)
    recs[0].screening_decision = "include"
    recs[0].screened_by = "AG Zippo"
    pipeline_02.write_screening_csv(recs, csv_path)

    pipeline_02.apply_auto_exclusions(csv_path, VOCAB, OMICS)

    a = pipeline_02.audit_sample(csv_path, n=5, seed=42)
    b = pipeline_02.audit_sample(csv_path, n=5, seed=42)
    assert [r.uid for r in a] == [r.uid for r in b]
    assert len(a) == 5
    assert all(r.screened_by == "auto" for r in a)


# --- appended: epigenomics/exome omics-term coverage (Task 8 review fix) ---
REAL_VOCAB = load_pain_vocabulary()


def test_atac_seq_of_drg_after_cci_is_left_for_a_human(pipeline_02):
    """Epigenomics family (ATAC-seq) was missing from DEFAULT_OMICS_TERMS,
    so a real chronic-pain study was wrongly auto-excluded as no_omics_term."""
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "ATAC-seq of dorsal root ganglion after chronic constriction injury"),
        REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)
    assert (decision, reason) == ("", "")


def test_genome_wide_methylation_in_fibromyalgia_is_left_for_a_human(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "genome-wide DNA methylation in fibromyalgia"),
        REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)
    assert (decision, reason) == ("", "")


def test_chip_seq_in_sni_model_is_left_for_a_human(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "ChIP-seq in the SNI model of neuropathic pain"),
        REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)
    assert (decision, reason) == ("", "")


def test_whole_exome_sequencing_survives_auto_exclusion(pipeline_02):
    """Bare 3-char acronyms "wes"/"wgs" are dropped from DEFAULT_OMICS_TERMS
    (they collide with common words, e.g. "wes" inside "western"). Coverage
    for exome/genome sequencing instead comes from "exome" (>=5 chars, safe
    as a substring) and the spelled-out "whole-exome"/"whole-genome" forms,
    so a real whole-exome study is still left for a human."""
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "whole-exome sequencing in migraine"),
        REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)
    assert (decision, reason) == ("", "")


def test_reset_auto_exclusions_clears_only_auto_rows(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    auto_excluded = _rec(pipeline_02, "maize leaf proteomics")
    auto_excluded.uid, auto_excluded.screening_decision = "a", "exclude"
    auto_excluded.exclusion_reason, auto_excluded.screened_by = "no_omics_term", "auto"

    human_included = _rec(pipeline_02, "RNA-seq of chronic pain in rats")
    human_included.uid, human_included.screening_decision = "b", "include"
    human_included.screened_by = "AG Zippo"

    human_excluded = _rec(pipeline_02, "acupuncture RCT for chronic pain")
    human_excluded.uid, human_excluded.screening_decision = "c", "exclude"
    human_excluded.exclusion_reason, human_excluded.screened_by = "wrong_condition", "AG Zippo"

    pipeline_02.write_screening_csv([auto_excluded, human_included, human_excluded], csv_path)

    n_reset = pipeline_02.reset_auto_exclusions(csv_path)
    assert n_reset == 1

    records = pipeline_02.load_existing_screening(csv_path)
    reset_rec = records["pubmed:a"]
    assert reset_rec.screening_decision == ""
    assert reset_rec.exclusion_reason == ""
    assert reset_rec.screened_by == ""

    kept_included = records["pubmed:b"]
    assert kept_included.screening_decision == "include"
    assert kept_included.screened_by == "AG Zippo"

    kept_excluded = records["pubmed:c"]
    assert kept_excluded.screening_decision == "exclude"
    assert kept_excluded.exclusion_reason == "wrong_condition"
    assert kept_excluded.screened_by == "AG Zippo"


def test_reset_preserves_human_override_with_stale_auto_tag(pipeline_02, tmp_path):
    """THE REGRESSION: a reviewer corrects an auto-exclusion to 'include' via
    --apply but only edits the decision column, leaving screened_by="auto"
    stale on the row. A later --reset-auto must NOT wipe that human decision."""
    csv_path = tmp_path / "screening.csv"
    stale_tag_include = _rec(pipeline_02, "RNA-seq of chronic pain in rats")
    stale_tag_include.uid = "stale"
    stale_tag_include.screening_decision = "include"
    stale_tag_include.screened_by = "auto"  # stale tag left by round-trip

    genuine_auto = _rec(pipeline_02, "maize leaf proteomics")
    genuine_auto.uid = "genuine"
    genuine_auto.screening_decision = "exclude"
    genuine_auto.exclusion_reason = "no_omics_term"
    genuine_auto.screened_by = "auto"

    pipeline_02.write_screening_csv([stale_tag_include, genuine_auto], csv_path)

    n_reset = pipeline_02.reset_auto_exclusions(csv_path)
    assert n_reset == 1

    records = pipeline_02.load_existing_screening(csv_path)
    preserved = records["pubmed:stale"]
    assert preserved.screening_decision == "include"
    assert preserved.screened_by == "auto"

    reset_rec = records["pubmed:genuine"]
    assert reset_rec.screening_decision == ""
    assert reset_rec.exclusion_reason == ""
    assert reset_rec.screened_by == ""


def test_reset_leaves_human_exclude_with_non_auto_reason_untouched(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    rec = _rec(pipeline_02, "acupuncture RCT for chronic pain")
    rec.uid = "h"
    rec.screening_decision = "exclude"
    rec.exclusion_reason = "wrong_species"
    rec.screened_by = "AG Zippo"
    pipeline_02.write_screening_csv([rec], csv_path)

    n_reset = pipeline_02.reset_auto_exclusions(csv_path)
    assert n_reset == 0

    result = pipeline_02.load_existing_screening(csv_path)["pubmed:h"]
    assert result.screening_decision == "exclude"
    assert result.exclusion_reason == "wrong_species"
    assert result.screened_by == "AG Zippo"


def test_reset_leaves_auto_tagged_row_with_non_auto_reason_untouched(pipeline_02, tmp_path):
    """A row bearing screened_by="auto" but a non-auto exclusion_reason (e.g.
    a human re-reasoned it as "duplicate" but left the stale tag) must be
    preserved, not reset."""
    csv_path = tmp_path / "screening.csv"
    rec = _rec(pipeline_02, "maize leaf proteomics")
    rec.uid = "reasoned"
    rec.screening_decision = "exclude"
    rec.exclusion_reason = "duplicate"
    rec.screened_by = "auto"
    pipeline_02.write_screening_csv([rec], csv_path)

    n_reset = pipeline_02.reset_auto_exclusions(csv_path)
    assert n_reset == 0

    result = pipeline_02.load_existing_screening(csv_path)["pubmed:reasoned"]
    assert result.screening_decision == "exclude"
    assert result.exclusion_reason == "duplicate"
    assert result.screened_by == "auto"


def test_cut_and_run_spelled_out_is_left_for_a_human(pipeline_02):
    decision, reason = pipeline_02.auto_screen(
        _rec(pipeline_02, "Cut and Run profiling of H3K27ac in DRG after CCI"),
        REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)
    assert (decision, reason) == ("", "")


def test_reset_then_reapply_rescues_previously_excluded_epigenomics_row(pipeline_02, tmp_path):
    csv_path = tmp_path / "screening.csv"
    rec = _rec(pipeline_02, "ATAC-seq of dorsal root ganglion after chronic constriction injury")
    rec.uid = "z"
    # Simulate the pre-fix state: auto-excluded under the old, narrower term list.
    rec.screening_decision = "exclude"
    rec.exclusion_reason = "no_omics_term"
    rec.screened_by = "auto"
    pipeline_02.write_screening_csv([rec], csv_path)

    n_reset = pipeline_02.reset_auto_exclusions(csv_path)
    assert n_reset == 1

    pipeline_02.apply_auto_exclusions(csv_path, REAL_VOCAB, pipeline_02.DEFAULT_OMICS_TERMS)

    result = pipeline_02.load_existing_screening(csv_path)["pubmed:z"]
    assert result.screening_decision == ""
    assert result.screened_by == ""
