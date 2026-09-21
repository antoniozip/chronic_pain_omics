"""Metabolomics Workbench discovery.

Every test here pins a defect that was live during development, not a
hypothetical: the single-hit response shape that reports one study as fifteen,
the truncated-stem regex that silently drops two thirds of the candidates, and
the sample-count sentinel that decides whether a study survives apply_filters.
"""
from __future__ import annotations

from cp_multiomics.ingest import metabolomics_workbench as wb

# A single-hit response is the study object itself, not a container.
SINGLE = {
    "study_id": "ST002974",
    "study_title": "Leishmania mexicana Promotes Pain-reducing Metabolomic Reprogramming",
    "species": "Mus musculus",
    "institute": "Ohio State University",
    "analysis_type": "LC-MS",
    "number_of_samples": "66",
    "submission_date": "2023-11-15",
    "release_date": "2024-10-18",
    "version": "1",
    "revision_no": "1",
    "revision_datetime": "-",
    "revision_comment": "-",
    "license": "CC BY 4.0",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
    "study_url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=PAIN",
}

MULTI = {
    "1": {"study_id": "ST001121", "study_title": "Urine metabolites in interstitial cystitis",
          "species": "Homo sapiens", "analysis_type": "LC-MS", "number_of_samples": "43",
          "submission_date": "2019-02-01", "institute": "X"},
    "2": {"study_id": "ST003177", "study_title": "Longitudinal Study in Rheumatoid Arthritis",
          "species": "Homo sapiens", "analysis_type": "LC-MS", "number_of_samples": "2492",
          "submission_date": "2024-06-01", "institute": "Y"},
}


class TestResponseShapes:
    def test_single_hit_is_one_study_not_fifteen_fields(self):
        """The bare object carries ~15 metadata fields. Counting values() would
        report 15 studies. This is the defect that made a first probe read
        `pain` as fifteen hits."""
        studies = wb._studies(SINGLE)
        assert len(studies) == 1
        assert studies[0]["study_id"] == "ST002974"

    def test_rank_keyed_object_yields_each_study(self):
        studies = wb._studies(MULTI)
        assert [s["study_id"] for s in studies] == ["ST001121", "ST003177"]

    def test_empty_payloads_yield_nothing(self):
        for empty in (None, {}, [], ""):
            assert wb._studies(empty) == []

    def test_entries_without_study_id_are_skipped(self):
        assert wb._studies({"1": {"note": "no id"}, "2": MULTI["1"]}) == [MULTI["1"]]


class TestPainScreening:
    def test_truncated_stems_match_full_words(self):
        """A trailing \\b after a truncated stem never matches: `\\barthrit\\b`
        excludes "arthritis". That bug cut a first screen from 23 to 7."""
        for title in ("Plasma Metabolome in Rheumatoid Arthritis",
                      "Effect of streptozotocin treatment on neuropathy",
                      "Metabolomic profiling in endometriosis patients",
                      "Nociceptive thresholds after nerve injury",
                      "Analgesic response in chronic low back pain"):
            assert wb.is_pain_related(title), title

    def test_unrelated_titles_are_rejected(self):
        for title in ("Elevation of ALDH1A1 Ameliorates Podocyte Injury",
                      "Rumen microbiota of Tibetan sheep",
                      "Metabolomics of hepatic steatosis"):
            assert not wb.is_pain_related(title), title

    def test_missing_title_is_not_pain_related(self):
        assert not wb.is_pain_related("")
        assert not wb.is_pain_related(None)


class TestSpeciesAndSamples:
    def test_accepted_species(self):
        for s in ("Homo sapiens", "  mus musculus ", "Rattus norvegicus"):
            assert wb.species_accepted(s), s

    def test_rejected_species(self):
        for s in ("Ovis aries", "Drosophila melanogaster", "", None):
            assert not wb.species_accepted(s), s

    def test_sample_count_parses_the_string_form(self):
        assert wb._n_samples("66") == 66
        assert wb._n_samples(" 2492 ") == 2492

    def test_unknown_sample_count_is_minus_one_not_zero(self):
        """0 reads as "fewer than min_samples" and drops the study; -1 passes
        through as unknown."""
        for raw in (None, "", "n/a", "0", 0):
            assert wb._n_samples(raw) == -1, raw


class TestRecord:
    def test_url_is_built_from_the_accession(self):
        """The payload's study_url echoes the query keyword
        (...StudyID=PAIN) instead of naming the study."""
        rec = wb.to_record(SINGLE, "metabolomics")
        assert "ST002974" in rec.url
        assert "PAIN" not in rec.url

    def test_record_fields(self):
        rec = wb.to_record(SINGLE, "metabolomics")
        assert rec.source == "METABOLOMICS_WORKBENCH"
        assert rec.accession == "ST002974"
        assert rec.species == ["Mus musculus"]
        assert rec.n_samples == 66
        assert rec.year == "2023"
        assert rec.modality == "metabolomics"

    def test_missing_species_gives_empty_list_not_empty_string(self):
        rec = wb.to_record({"study_id": "ST1", "study_title": "t"}, "metabolomics")
        assert rec.species == []
        assert rec.n_samples == -1


class TestScreen:
    def test_every_excluded_study_carries_a_reason_code(self):
        """PRISMA: never silently drop a study."""
        studies = [
            SINGLE,                                                    # keep
            {"study_id": "ST9", "study_title": "Hepatic steatosis",
             "species": "Homo sapiens"},                               # no pain
            {"study_id": "ST8", "study_title": "Arthritis in sheep",
             "species": "Ovis aries"},                                 # wrong species
        ]
        kept, excluded = wb.screen(studies)
        assert [s["study_id"] for s in kept] == ["ST002974"]
        assert {e["accession"]: e["reason_code"] for e in excluded} == {
            "ST9": "no_pain_phenotype", "ST8": "wrong_species"}
        assert len(kept) + len(excluded) == len(studies)

    def test_lipidomics_arm_inverts_the_modality_test(self):
        lipid = {"study_id": "ST001273", "species": "Mus musculus",
                 "study_title": "Lipidomics Dataset of Traumatic Optic Neuropathy in Mice"}
        kept_m, excluded_m = wb.screen([lipid], modality="metabolomics")
        kept_l, _ = wb.screen([lipid], modality="lipidomics")
        assert kept_m == []
        assert excluded_m[0]["reason_code"] == "wrong_modality_arm"
        assert [s["study_id"] for s in kept_l] == ["ST001273"]

    def test_a_study_is_claimed_by_exactly_one_arm(self):
        studies = [SINGLE, {"study_id": "ST001273", "species": "Mus musculus",
                            "study_title": "Lipidomics of Optic Neuropathy"}]
        m, _ = wb.screen(studies, modality="metabolomics")
        lp, _ = wb.screen(studies, modality="lipidomics")
        assert set(s["study_id"] for s in m) & set(s["study_id"] for s in lp) == set()
