from cp_multiomics.genomics.resolve import (
    SumstatsStudy,
    normalize_phenotype,
    resolve_sumstats_studies,
)

_STUDIES = {
    "back pain": [
        {"accessionId": "GCST001", "fullPvalueSet": True,
         "diseaseTrait": {"trait": "Back pain"},
         "initialSampleSize": "5,000 UK Biobank cases, 100,000 controls"},
        {"accessionId": "GCST002", "fullPvalueSet": False,
         "diseaseTrait": {"trait": "Back pain (no sumstats)"},
         "initialSampleSize": "10 cases"},
    ],
    "fibromyalgia": [
        {"accessionId": "GCST003", "fullPvalueSet": True,
         "diseaseTrait": {"trait": "Fibromyalgia"},
         "initialSampleSize": "1,362 European cases, 44,047 controls"},
        # same accession also appears under back pain query -> dedup
        # (trait deliberately differs from the "back pain" entry so
        # first-wins vs last-wins is observable)
        {"accessionId": "GCST001", "fullPvalueSet": True,
         "diseaseTrait": {"trait": "Back pain (secondary annotation)"},
         "initialSampleSize": "5,000 UK Biobank cases, 100,000 controls"},
    ],
}


def _fake_get(url, params=None):
    trait = params["efoTrait"]
    return {"_embedded": {"studies": _STUDIES.get(trait, [])}}


def test_normalize_phenotype_strips_parens_and_case():
    assert normalize_phenotype("Knee pain for three months (UKB data field 3773)") \
        == "knee pain for three months"
    assert normalize_phenotype("Back pain") == "back pain"


def test_only_full_pvalue_set_studies_are_kept():
    out = resolve_sumstats_studies(["back pain"], _fake_get)
    accs = {s.accession for s in out}
    assert accs == {"GCST001"}          # GCST002 dropped (fullPvalueSet False)


def test_accession_deduplicated_across_traits():
    out = resolve_sumstats_studies(["back pain", "fibromyalgia"], _fake_get)
    accs = [s.accession for s in out]
    assert sorted(accs) == ["GCST001", "GCST003"]
    assert len(accs) == len(set(accs))  # no duplicate GCST001


def test_study_fields_are_parsed():
    out = {s.accession: s for s in resolve_sumstats_studies(["fibromyalgia"], _fake_get)}
    s = out["GCST003"]
    assert s.n == 45409
    assert s.cohort == ""
    assert s.is_burden is False
    assert s.phenotype == "fibromyalgia"
    assert isinstance(s, SumstatsStudy)


def test_ukb_cohort_detected_from_sample_size():
    out = {s.accession: s for s in resolve_sumstats_studies(["back pain"], _fake_get)}
    assert out["GCST001"].cohort == "UKB"


def test_normalize_phenotype_falls_back_when_all_parenthetical():
    assert normalize_phenotype("(UKB data field 3773)") == "ukb data field 3773"


def test_normalize_phenotype_distinguishes_two_all_parenthetical_traits():
    # two different annotation-only traits must NOT collapse to the same key
    a = normalize_phenotype("(UKB data field 3773)")
    b = normalize_phenotype("(UKB data field 3414)")
    assert a != b
    assert a and b   # neither is empty


def test_dedup_keeps_first_occurrence_content():
    # GCST001 appears under "back pain" (queried first) and "fibromyalgia";
    # the retained record must carry the FIRST trait, not the later duplicate.
    out = {
        s.accession: s
        for s in resolve_sumstats_studies(["back pain", "fibromyalgia"], _fake_get)
    }
    assert out["GCST001"].trait == "Back pain"
