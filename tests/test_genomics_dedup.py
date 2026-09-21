from cp_multiomics.genomics.dedup import (
    PRESUMED_BIOBANK,
    _canonical_cohort,
    apply_presumed_cohort,
    dedup_studies,
)
from cp_multiomics.genomics.resolve import SumstatsStudy


def _s(acc, n, cohort, phenotype, burden=False):
    return SumstatsStudy(accession=acc, trait=phenotype, sample_size_text="",
                         n=n, cohort=cohort, is_burden=burden, phenotype=phenotype)


def test_burden_studies_are_dropped_with_reason():
    kept, dropped = dedup_studies([
        _s("GCST1", 1000, "", "knee pain"),
        _s("GCST2", 2000, "UKB", "knee pain", burden=True),
    ])
    assert [s.accession for s in kept] == ["GCST1"]
    assert ("GCST2", "burden") in [(s.accession, r) for s, r in dropped]


def test_same_cohort_and_phenotype_keeps_largest_n():
    kept, dropped = dedup_studies([
        _s("GCST1", 50000, "UKB", "back pain"),
        _s("GCST2", 90000, "UKB", "back pain"),
        _s("GCST3", 10000, "UKB", "back pain"),
    ])
    assert [s.accession for s in kept] == ["GCST2"]           # largest N wins
    reasons = {(s.accession, r) for s, r in dropped}
    assert reasons == {("GCST1", "cohort_overlap"), ("GCST3", "cohort_overlap")}


def test_same_cohort_different_phenotype_both_kept():
    kept, _ = dedup_studies([
        _s("GCST1", 50000, "UKB", "knee pain"),
        _s("GCST2", 50000, "UKB", "hip pain"),
    ])
    assert sorted(s.accession for s in kept) == ["GCST1", "GCST2"]


def test_empty_cohort_studies_are_never_collapsed():
    # two no-cohort studies on the same phenotype are presumed independent
    kept, dropped = dedup_studies([
        _s("GCST1", 500, "", "fibromyalgia"),
        _s("GCST2", 800, "", "fibromyalgia"),
    ])
    assert sorted(s.accession for s in kept) == ["GCST1", "GCST2"]
    assert dropped == []


def test_tie_on_n_breaks_by_accession():
    kept, _ = dedup_studies([
        _s("GCSTB", 1000, "UKB", "back pain"),
        _s("GCSTA", 1000, "UKB", "back pain"),
    ])
    assert [s.accession for s in kept] == ["GCSTA"]           # lexicographic tiebreak


def test_burden_study_within_overlap_group_does_not_pollute_cohort_dedup():
    kept, dropped = dedup_studies([
        _s("GCSTB", 9000, "UKB", "back pain", burden=True),
        _s("GCST1", 3000, "UKB", "back pain"),
        _s("GCST2", 4000, "UKB", "back pain"),
    ])
    assert [s.accession for s in kept] == ["GCST2"]
    dropped_dict = dict((s.accession, r) for s, r in dropped)
    assert dropped_dict == {"GCSTB": "burden", "GCST1": "cohort_overlap"}


def test_tiebreak_smallest_accession_with_differing_lengths():
    # "GCST1" must beat "GCST10" on equal n (lexicographic min), which the old
    # negated-ord list trick got wrong.
    kept, _ = dedup_studies([
        _s("GCST10", 1000, "UKB", "knee pain"),
        _s("GCST1", 1000, "UKB", "knee pain"),
    ])
    assert [s.accession for s in kept] == ["GCST1"]


def test_count_preservation_every_study_kept_or_dropped_once():
    studies = [
        _s("GCSTB", 9000, "UKB", "back pain", burden=True),
        _s("GCST1", 3000, "UKB", "back pain"),
        _s("GCST2", 4000, "UKB", "back pain"),
        _s("GCST3", 500, "", "fibromyalgia"),
    ]
    kept, dropped = dedup_studies(studies)
    assert len(kept) + len(dropped) == len(studies)
    all_accs = {s.accession for s in kept} | {s.accession for s, _ in dropped}
    assert all_accs == {s.accession for s in studies}  # none lost or duplicated


def test_large_n_empty_cohort_becomes_presumed_biobank():
    studies = [_s("GCST1", 400000, "", "back pain"), _s("GCST2", 500, "", "fibromyalgia")]
    out = {s.accession: s for s in apply_presumed_cohort(studies)}
    assert out["GCST1"].cohort == PRESUMED_BIOBANK      # N>100k -> presumed
    assert out["GCST2"].cohort == ""                    # small N, stays independent


def test_icd10_phenotype_empty_cohort_becomes_presumed_biobank():
    out = {s.accession: s for s in apply_presumed_cohort(
        [_s("GCST1", 5000, "", "icd10 m54: dorsalgia")])}   # small N but ICD10-coded
    assert out["GCST1"].cohort == PRESUMED_BIOBANK


def test_real_cohort_is_not_overwritten():
    out = {s.accession: s for s in apply_presumed_cohort(
        [_s("GCST1", 400000, "UKB", "back pain")])}
    assert out["GCST1"].cohort == "UKB"                 # already detected, unchanged


def test_presumed_biobank_then_dedup_collapses_duplicate_phenotypes():
    # end-to-end: 3 large-N empty-cohort back-pain studies collapse to 1
    studies = apply_presumed_cohort([
        _s("GCST1", 200000, "", "back pain"),
        _s("GCST2", 400000, "", "back pain"),
        _s("GCST3", 300000, "", "back pain"),
        _s("GCST4", 350000, "", "fibromyalgia"),   # different phenotype -> kept
    ])
    kept, dropped = dedup_studies(studies)
    kept_accs = sorted(s.accession for s in kept)
    assert kept_accs == ["GCST2", "GCST4"]          # max-N back pain + the fibromyalgia
    assert len(dropped) == 2                          # GCST1, GCST3 as cohort_overlap


def test_canonical_cohort_maps_presumed_to_ukb():
    assert _canonical_cohort(PRESUMED_BIOBANK) == "UKB"
    assert _canonical_cohort("UKB") == "UKB"
    assert _canonical_cohort("FinnGen") == "FinnGen"
    assert _canonical_cohort("") == ""


def test_ukb_and_presumed_biobank_same_phenotype_collapse_together():
    # one literally-UKB and one presumed-biobank back-pain study must collapse
    kept, dropped = dedup_studies([
        _s("GCST1", 300000, "UKB", "back pain"),
        _s("GCST2", 500000, PRESUMED_BIOBANK, "back pain"),
    ])
    assert [s.accession for s in kept] == ["GCST2"]          # max N wins across the merged bucket
    assert dropped[0][0].accession == "GCST1" and dropped[0][1] == "cohort_overlap"


def test_kept_study_retains_original_cohort_label_for_audit():
    # the surviving study keeps PRESUMED_BIOBANK (not rewritten to UKB) so the
    # named-vs-inferred split stays visible in the manifest
    kept, _ = dedup_studies([
        _s("GCST1", 300000, "UKB", "back pain"),
        _s("GCST2", 500000, PRESUMED_BIOBANK, "back pain"),
    ])
    assert kept[0].cohort == PRESUMED_BIOBANK


def test_finngen_not_merged_into_ukb():
    # a FinnGen study on the same phenotype must NOT collapse with a UKB one
    kept, _ = dedup_studies([
        _s("GCST1", 300000, "UKB", "back pain"),
        _s("GCST2", 300000, "FinnGen", "back pain"),
    ])
    assert sorted(s.accession for s in kept) == ["GCST1", "GCST2"]
