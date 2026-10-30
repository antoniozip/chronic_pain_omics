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


# --- appended: B7 -- a failed streamed download must not leak its temp file --

def test_failed_streamed_download_leaves_no_temp_file(tmp_path, monkeypatch):
    class _Head:
        headers = {"Content-Length": str(500 * 1024 * 1024)}   # forces disk path

    def _boom(*args, **kwargs):
        raise OSError("404 Not Found")

    monkeypatch.setattr(pipeline_05d.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(pipeline_05d.requests, "head", lambda *a, **k: _Head())
    monkeypatch.setattr(pipeline_05d.requests, "get", _boom)

    assert pipeline_05d.fetch_tarball("https://example/GSE1_RAW.tar") is None
    assert list(tmp_path.glob("geo_tar_*")) == []



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


# -- pseudo-aligner outputs ----------------------------------------------------
# GSE221957 (Salmon) and GSE303635 (Salmon, transcript level) were merged with
# every sample holding the same value for every gene: the generic branch of
# parse_sample_content takes the first numeric column, which in a Salmon
# quant.sf is `Length`, the same in every sample. Step 05 then stopped with
# "all genes have equal values". Both deposits sat in the sheet as an
# unexplained step-05 error. The wrong column is silent unless the matrix is
# inspected, so the guard is the choice of column by header.

SALMON = (
    "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
    "ENSRNOT00000094517.1\t303\t131.000\t0.000000\t0.000\n"
    "ENSRNOT00000064892.3\t402\t253.159\t1.642825\t8.764\n"
    "ENSRNOT00000095068.1\t144\t24.000\t30.188045\t15.267\n"
)
KALLISTO = (
    "target_id\tlength\teff_length\test_counts\ttpm\n"
    "ENST00000641515\t2618\t2400.0\t12.5\t0.4\n"
    "ENST00000426406\t939\t720.0\t301.0\t9.9\n"
)


def test_a_salmon_quant_file_yields_numreads_not_the_transcript_length() -> None:
    s = pipeline_05d.parse_sample_content(SALMON)
    assert s["ENSRNOT00000064892.3"] == 8.764
    assert s["ENSRNOT00000095068.1"] == 15.267
    assert s["ENSRNOT00000094517.1"] == 0.0


def test_a_kallisto_abundance_file_yields_est_counts_not_the_length() -> None:
    s = pipeline_05d.parse_sample_content(KALLISTO)
    assert s["ENST00000641515"] == 12.5
    assert s["ENST00000426406"] == 301.0


def test_two_salmon_samples_do_not_come_out_identical() -> None:
    other = SALMON.replace("8.764", "40.100").replace("15.267", "2.000")
    a = pipeline_05d.parse_sample_content(SALMON)
    b = pipeline_05d.parse_sample_content(other)
    assert not a.equals(b)


def test_a_plain_two_column_count_file_is_unaffected() -> None:
    s = pipeline_05d.parse_sample_content("gene\tcount\nA\t5\nB\t0\n")
    assert s["A"] == 5.0 and s["B"] == 0.0
