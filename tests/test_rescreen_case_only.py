"""Cover the control-detection the re-screen turns on.

The automatic screen missed a body of endometriosis case/control studies
because `conf/analysis/default.yaml` names only pain terms. Every check here
guards the two ways the replacement can fail the same silent way: calling a
case a control, or ranking a cycle-phase field above a disease field and
proposing the wrong contrast.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import rescreen_case_only_studies as rs  # noqa: E402


@pytest.mark.parametrize("value", [
    "control", "healthy", "normal", "healthy control", "no endometriosis",
    "non-endometriosis", "non endometriosis", "disease free endometrium",
    "unaffected", "sham", "naive", "vehicle",
])
def test_control_like_values_are_recognised(value: str):
    assert rs.is_control(value)


@pytest.mark.parametrize("value", [
    "endometriosis", "endometriosis stage iv", "endometriotic tissue",
    "interstitial cystitis", "fibromyalgia", "",
])
def test_case_like_values_are_not_mistaken_for_controls(value: str):
    # "non-endometriosis" contains "endometriosis", so the two directions have
    # to be tested separately: a rule loose enough to catch the control label
    # will catch the case label too unless it is anchored.
    assert not rs.is_control(value)


def test_split_counts_each_side_and_names_it():
    values = ["endometriosis"] * 3 + ["no endometriosis"] * 2
    n_ctrl, n_case, case_label, ctrl_label = rs.split_on(values)
    assert (n_ctrl, n_case) == (2, 3)
    assert case_label == "endometriosis"
    assert ctrl_label == "no endometriosis"


def test_a_field_with_no_control_value_is_not_a_candidate():
    chars = {"disease": ["endometriosis"] * 4}
    assert rs.candidate_fields(chars, 4) == []


def test_a_field_with_no_case_value_is_not_a_candidate():
    chars = {"disease": ["control"] * 4}
    assert rs.candidate_fields(chars, 4) == []


def test_a_field_whose_length_does_not_match_the_samples_is_skipped():
    # A SuperSeries part can contribute a field to only some of its samples.
    chars = {"disease": ["endometriosis", "control"]}
    assert rs.candidate_fields(chars, 4) == []


def test_the_disease_field_outranks_a_cycle_phase_field():
    chars = {
        "menstrual cycle phase": ["proliferative", "control", "proliferative",
                                  "control"],
        "disease": ["endometriosis", "no endometriosis", "endometriosis",
                    "no endometriosis"],
    }
    ranked = rs.candidate_fields(chars, 4)
    assert ranked[0]["field"] == "disease"
    assert ranked[0]["looks_non_disease"] is False


@pytest.mark.parametrize("field", [
    "menstrual cycle phase", "tissue", "cell type", "treatment", "age", "sex",
    "batch", "timepoint", "donor",
])
def test_fields_that_separate_without_describing_disease_are_flagged(field: str):
    chars = {field: ["treated", "control", "treated", "control"]}
    found = rs.candidate_fields(chars, 4)
    assert found and found[0]["looks_non_disease"] is True


def test_the_audit_reports_rather_than_decides():
    # A separating field may describe cycle phase or a knockdown rather than
    # disease, which is a reading of the study. The module must not expose a
    # verdict, or a reviewer will take one.
    assert not hasattr(rs, "verdict")
    assert "n_case" in rs.COLUMNS and "verdict" not in rs.COLUMNS


# -- negation binds to a diagnosis, not to a feature of one ------------------
# GSE28242 deposits three arms: normal, "PBS without lesions" and "PBS with
# lesions". A bare `without` read the middle arm as control and put five
# painful bladder syndrome patients in the control group, which is a
# well-formed split and the wrong one.

@pytest.mark.parametrize("value", [
    "pbs without lesions", "with lesions", "lesion-free", "without stimulation",
    "no stimulus",
])
def test_negating_a_feature_does_not_make_a_control(value: str):
    assert not rs.is_control(value)


@pytest.mark.parametrize("value", [
    "without endometriosis", "no endometriosis", "without disease",
    "no pain", "non-neuropathic", "disease-free", "pain-free",
])
def test_negating_a_diagnosis_does_make_a_control(value: str):
    assert rs.is_control(value)


def test_the_gse28242_arms_split_the_way_the_deposit_states_them():
    values = (["normal"] * 5 + ["pbs without lesions"] * 5
              + ["pbs with lesions"] * 3)
    n_ctrl, n_case, _, _ = rs.split_on(values)
    assert (n_ctrl, n_case) == (5, 8)


@pytest.mark.parametrize("field", [
    "disease stage", "disease status", "endometriosis stage", "diagnosis",
])
def test_a_disease_field_is_not_flagged_by_a_substring(field: str):
    # Unanchored `age` matches inside "stage". That flagged GSE141549's
    # "disease stage" as non-disease and let a cycle-phase field outrank it on
    # the largest study in the target list, 408 samples.
    chars = {field: ["endometriosis", "no disease", "endometriosis", "no disease"]}
    found = rs.candidate_fields(chars, 4)
    assert found and found[0]["looks_non_disease"] is False


# -- the curated rules for studies the screen cannot do ----------------------

import assign_rescued_groups as arg  # noqa: E402


def _needs_series_matrix(accession: str) -> None:
    """Skip when the deposit's series matrix is not in the local GEO cache.

    `arg.assign` reads `data/raw/geo_cache/`, which is gitignored, so on a
    fresh clone and in CI it finds nothing and returns an empty assignment.
    Without this the tests below do not skip, they *fail*, asserting 0 == 12
    against an absent input -- which says nothing about the curated rule and
    turns every CI run red on a checkout that never had the data.
    """
    if not (list(arg.CACHE.glob(f"{accession}-*_series_matrix.txt.gz"))
            + list(arg.CACHE.glob(f"{accession}_series_matrix.txt.gz"))):
        pytest.skip(f"no series matrix cached for {accession}; "
                    "run scripts/download_series_matrices.py")


def test_every_curated_rule_carries_a_note():
    # Reading these rules is how the assignment gets reviewed; there is no
    # other record of why a given arm was chosen.
    for accession, rule in arg.CURATED.items():
        assert rule.get("note"), accession
        assert len(rule["note"]) > 40, accession


def test_every_curated_rule_names_both_arms():
    for accession, rule in arg.CURATED.items():
        has_lists = "case" in rule and "control" in rule
        has_patterns = "case_match" in rule and "control_match" in rule
        assert has_lists or has_patterns, accession


def test_a_curated_rule_drops_what_it_does_not_name():
    # GSE621 deposits two designs; the 16 antiproliferative-factor samples must
    # not enter as cases. GSE11783's five non-ulcer biopsies come from the same
    # patients as its ulcer ones. GSE269047's ME/CFS arm is not a pain arm.
    for accession, expected in (("GSE621", 12), ("GSE11783", 11),
                                ("GSE269047", 35)):
        _needs_series_matrix(accession)
        assert len(arg.assign(accession, arg.CURATED[accession])) == expected


def test_no_curated_study_is_also_in_the_single_cell_arm():
    import csv as _csv
    with open(REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv") as fh:
        single_cell = {r["accession"] for r in _csv.DictReader(fh)}
    assert not (set(arg.CURATED) & single_cell)


def test_gse51981_excludes_controls_with_their_own_pelvic_pathology():
    # The automatic screen splits it 77 v 71 on `endometriosis/no
    # endometriosis`, and 37 of those 71 controls carry uterine pelvic
    # pathology -- fibroids or adenomyosis, themselves causes of pelvic pain.
    # Over half the control arm would be a pain population.
    _needs_series_matrix("GSE51981")
    pairs = arg.assign("GSE51981", arg.CURATED["GSE51981"])
    n_case = sum(1 for _, g in pairs if g == "case")
    assert (n_case, len(pairs) - n_case) == (77, 34)
    assert len(pairs) == 111
