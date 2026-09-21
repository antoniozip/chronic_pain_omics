"""Cover the two defects that hid case/control splits, and the guards around them.

`auto_annotate_geo` set every one of the 75 `no_case_control_split_in_series_matrix`
verdicts in the human arm. Two defects put them there, and both produced the
same silent output: a study with a perfectly good control arm reporting no
split at all, indistinguishable from a study that has none.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "pipeline"))

import auto_annotate_geo as aag  # noqa: E402

# -- defect 1: a control that names the disease it does not have -------------

@pytest.mark.parametrize("value", [
    "no endometriosis", "non-endometriosis", "non endometriosis",
    "non_endometriosis", "without endometriosis", "endometriosis-free",
    "disease free endometrium", "disease-free", "no disease",
])
def test_a_negated_disease_is_a_control(value: str):
    # Before the fix these classified as case, because the disease token is
    # inside its own negation. Both arms became cases and the split vanished.
    assert aag._classify_value(value) == "control"


@pytest.mark.parametrize("value", [
    "endometriosis", "endometriosis stage iv", "neuropathic pain",
    "interstitial cystitis", "fibromyalgia",
])
def test_the_disease_itself_is_still_a_case(value: str):
    assert aag._classify_value(value) == "case"


def test_the_gse51981_arms_split():
    values = ["endometriosis"] * 3 + ["non-endometriosis"] * 2
    labels = [aag._classify_value(v) for v in values]
    assert labels.count("case") == 3
    assert labels.count("control") == 2


# -- defect 2: a skip token matching inside a longer word --------------------

@pytest.mark.parametrize("key", [
    "disease stage", "disease", "disease status", "disease state",
    "endometriosis status", "sample group", "disease severity", "diagnosis",
])
def test_a_disease_column_is_not_skipped(key: str):
    # "age" is inside "stage". That discarded the disease column of GSE141549,
    # the largest study in the arm.
    assert not aag._column_should_skip(key)


@pytest.mark.parametrize("key", [
    "age", "sex", "tissue", "cell type", "batch", "donor", "platform",
])
def test_a_technical_column_is_still_skipped(key: str):
    assert aag._column_should_skip(key)


@pytest.mark.parametrize("key", [
    "biopsied site", "biopsy site", "anatomical location", "region",
])
def test_an_anatomical_column_is_skipped(key: str):
    # GSE238208 biopsies Hunner lesions and non-lesional mucosa from the same
    # bladders. Read as a contrast, a within-patient site comparison becomes
    # disease against control.
    assert aag._column_should_skip(key)


def test_a_lesion_type_is_not_a_diagnosis():
    # "hunner" was added to the vocabulary and removed again: with the negation
    # rule it manufactured a 25 v 25 contrast out of lesional and non-lesional
    # biopsies of the same disease.
    assert "hunner" not in [t.lower() for t in aag.PAIN_TOKENS]


# -- the conditions the search knows about -----------------------------------

@pytest.mark.parametrize("condition", [
    "endometriosis", "interstitial cystitis", "low back pain", "sciatica",
    "temporomandibular", "fibromyalgia", "osteoarthritis",
])
def test_the_amendment_conditions_are_in_the_vocabulary(condition: str):
    # The 2026-08-28 amendment added these to the search query. A condition the
    # retrieval can find and the annotation cannot is a study retrieved and
    # then dropped for having no split.
    assert condition in [t.lower() for t in aag.PAIN_TOKENS]


# -- the single-cell guard ---------------------------------------------------

def test_the_single_cell_arm_is_reserved():
    # Every single-cell accession is also in the transcriptomic screening
    # record, and {accession}_groups.csv in the shared cache is what
    # run_transcriptomics_da() reads. One written here wires the arms together.
    reserved = aag.single_cell_accessions()
    assert reserved
    with open(REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv") as fh:
        expected = {row["accession"] for row in csv.DictReader(fh)}
    assert reserved == expected


def test_a_reserved_accession_gets_no_groups_file(tmp_path: Path):
    accession = sorted(aag.single_cell_accessions())[0]
    (tmp_path / f"{accession}_series_matrix.txt.gz").write_bytes(b"")
    aag.run(tmp_path, dry_run=False)
    assert not (tmp_path / f"{accession}_groups.csv").exists()


# -- the dry run -------------------------------------------------------------

def test_a_dry_run_writes_nothing_into_the_cache(tmp_path: Path):
    # Re-running the screen over the real cache is how an audit turns into a
    # change. The dry run has to leave the tree exactly as it found it.
    (tmp_path / "GSE1_series_matrix.txt.gz").write_bytes(b"")
    before = sorted(p.name for p in tmp_path.iterdir())
    aag.run(tmp_path, dry_run=True)
    assert sorted(p.name for p in tmp_path.iterdir()) == before
