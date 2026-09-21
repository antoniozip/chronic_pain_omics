import pandas as pd
import pytest

from cp_multiomics.genomics.grouping import (
    assign_pain_group,
    poolable_groups,
    select_representatives,
)


@pytest.mark.parametrize("phenotype,expected", [
    ("icd10 m79.7: fibromyalgia", "fibromyalgia"),
    ("Fibromyalgia", "fibromyalgia"),
    ("chronic back pain", "back_pain"),
    ("icd10 m54.4: lumbago with sciatica", "back_pain"),
    ("chronic knee pain", "knee_pain"),               # knee before back
    ("chronic hip pain", "hip_pain"),                 # hip before back
    ("central neuropathic pain", "neuropathic_pain"),
    ("postlaminectomy syndrome", "neuropathic_pain"),
    ("chronic widespread pain", "chronic_widespread_pain"),  # before fibromyalgia
    ("degree bothered by headaches in the last three months - 602 bothered a lot", "headache"),
    ("icd10 r52.2: other chronic pain", "chronic_pain_broad"),
    ("icd10 r51: headache", "headache"),
    ("menstruation quality of life impact", "nonpain"),
    ("response to bnt162b2 vaccine", "nonpain"),
    ("headache after covid-19 booster vaccination", "nonpain"),  # adverse event beats headache
])
def test_assign_pain_group(phenotype, expected):
    assert assign_pain_group(phenotype) == expected


def test_unmapped_phenotype():
    assert assign_pain_group("some unrelated trait") == "unmapped"


def _df(rows):
    return pd.DataFrame(rows)


def test_ukb_slices_collapse_to_max_n_representative():
    df = _df([
        {"accession": "GCST_A", "n": 100000, "cohort": "UKB", "phenotype": "back pain"},
        {"accession": "GCST_B", "n": 200000, "cohort": "PRESUMED_BIOBANK",
         "phenotype": "icd10 m54: dorsalgia"},
    ])
    out = select_representatives(df).set_index("accession")
    assert out.loc["GCST_B", "is_representative"]           # larger n wins
    assert out.loc["GCST_B", "pooled_status"] == "representative"
    assert not out.loc["GCST_A", "is_representative"]
    assert out.loc["GCST_A", "pooled_status"] == "collapsed_overlap"
    # both are the same UKB cohort family
    assert out.loc["GCST_A", "cohort_family"] == out.loc["GCST_B", "cohort_family"] == "UKB"


def test_empty_cohort_studies_are_each_independent():
    df = _df([
        {"accession": "GCST_C", "n": 5000, "cohort": "", "phenotype": "fibromyalgia"},
        {"accession": "GCST_D", "n": 6000, "cohort": None, "phenotype": "icd10 m79: fibromyalgia"},
    ])
    out = select_representatives(df).set_index("accession")
    assert out.loc["GCST_C", "is_representative"]
    assert out.loc["GCST_D", "is_representative"]
    assert out.loc["GCST_C", "cohort_family"] != out.loc["GCST_D", "cohort_family"]


def test_nonpain_never_representative():
    df = _df([
        {"accession": "GCST_E", "n": 5734, "cohort": "",
         "phenotype": "menstruation quality of life impact"},
    ])
    out = select_representatives(df).set_index("accession")
    assert out.loc["GCST_E", "pooled_status"] == "nonpain"
    assert not out.loc["GCST_E", "is_representative"]


def test_poolable_groups_needs_two_cohorts():
    df = _df([
        {"accession": "A", "n": 100000, "cohort": "UKB", "phenotype": "back pain"},
        {"accession": "B", "n": 50000, "cohort": "", "phenotype": "chronic back pain"},
        {"accession": "C", "n": 200000, "cohort": "UKB", "phenotype": "headaches"},
    ])
    rep = select_representatives(df)
    counts = poolable_groups(rep)
    assert counts["back_pain"] == 2      # UKB + one independent
    assert counts["headache"] == 1       # single UKB cohort
    poolable = {g: c for g, c in counts.items() if c >= 2}
    assert "back_pain" in poolable
    assert "headache" not in poolable
