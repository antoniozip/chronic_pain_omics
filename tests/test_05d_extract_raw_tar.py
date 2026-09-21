"""Tests for the _RAW.tar extraction step.

The defect these cover: GSE205494's tarball carries both an FPKM file and a
raw-count file for every GSM, and the merge loop added a column per *file*.
The result was a 24-column matrix from 12 samples whose columns alternated
between two quantification scales -- 0.27, 9.0, 0.12, 4.0 -- which then fails
`is_integer_counts` and goes through limma-trend as if it were one homogeneous
assay, with the sample count doubled. See plan/2026-08-29-amendment-da-run.md.
"""

import importlib.util
import tarfile
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_05d", Path(__file__).parent.parent / "pipeline" / "05d_extract_raw_tar.py"
)
pipeline_05d = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_05d)


def _members(*names: str) -> list[tarfile.TarInfo]:
    out = []
    for n in names:
        info = tarfile.TarInfo(name=n)
        info.size = 10
        out.append(info)
    return out


def test_one_member_per_gsm_prefers_raw_counts() -> None:
    members = _members(
        "GSM1_sciatic_FPKM.txt.gz",
        "GSM1_sciatic_raw_counts.txt.gz",
        "GSM2_sciatic_FPKM.txt.gz",
        "GSM2_sciatic_raw_counts.txt.gz",
    )
    kept = pipeline_05d.select_one_member_per_sample(members)
    assert [m.name for m in kept] == [
        "GSM1_sciatic_raw_counts.txt.gz",
        "GSM2_sciatic_raw_counts.txt.gz",
    ]


def test_a_sample_with_one_file_is_untouched() -> None:
    members = _members("GSM1_counts.txt.gz", "GSM2_counts.txt.gz")
    kept = pipeline_05d.select_one_member_per_sample(members)
    assert [m.name for m in kept] == ["GSM1_counts.txt.gz", "GSM2_counts.txt.gz"]


def test_a_sample_with_only_a_derived_quantification_is_kept() -> None:
    # Declining here would silently drop the sample; FPKM still analyses, via
    # limma-trend rather than DESeq2.
    members = _members("GSM1_TPM.txt.gz", "GSM2_TPM.txt.gz")
    kept = pipeline_05d.select_one_member_per_sample(members)
    assert len(kept) == 2


def test_selection_is_deterministic_when_nothing_distinguishes_the_files() -> None:
    members = _members("GSM1_b.txt.gz", "GSM1_a.txt.gz")
    first = pipeline_05d.select_one_member_per_sample(members)
    second = pipeline_05d.select_one_member_per_sample(list(reversed(members)))
    assert [m.name for m in first] == [m.name for m in second]
    assert len(first) == 1


def test_members_without_a_gsm_are_passed_through() -> None:
    # The auto-detection fallback downstream keys on filenames, so a tarball
    # whose members are not GSM-named must not be emptied here.
    members = _members("sampleA_counts.txt.gz", "sampleB_counts.txt.gz")
    kept = pipeline_05d.select_one_member_per_sample(members)
    assert len(kept) == 2


def test_a_name_saying_both_counts_and_fpkm_is_treated_as_derived() -> None:
    members = _members("GSM1_FPKM_counts.txt.gz", "GSM1_htseq.txt.gz")
    kept = pipeline_05d.select_one_member_per_sample(members)
    assert [m.name for m in kept] == ["GSM1_htseq.txt.gz"]


def test_the_streamed_tarball_survives_until_the_fallbacks_are_done(tmp_path) -> None:
    # The groups-based pass used to delete the temp file immediately after its
    # own loop, and both auto-detection fallbacks then re-opened a path that no
    # longer existed. Every study whose groups.csv GSM ids did not match its
    # tarball died on FileNotFoundError against its own download.
    src = Path(__file__).parent.parent / "pipeline" / "05d_extract_raw_tar.py"
    body = src.read_text()
    fallback_at = body.index("trying auto-detection")
    close_at = body.rindex("    tf.close()", 0, fallback_at)
    between = body[close_at:fallback_at]
    assert "os.unlink" not in between

    # and the helper actually removes the file when called
    doomed = tmp_path / "geo_tar_x.tar"
    doomed.write_bytes(b"x")
    pipeline_05d._discard_tarball(str(doomed))
    assert not doomed.exists()
    pipeline_05d._discard_tarball(str(doomed))  # idempotent
    pipeline_05d._discard_tarball(None)
