"""Tests for pipeline/05c_download_rnaseq_counts.py (B3).

An interrupted download used to stream straight into the final path, and
`already_downloaded()` only asked whether *some* count-patterned file existed,
so a truncated matrix was treated as a completed study forever.
"""

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_05c",
    Path(__file__).parent.parent / "pipeline" / "05c_download_rnaseq_counts.py",
)
pipeline_05c = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_05c)


class _Response:
    def __init__(self, headers, chunks, fail_after=None):
        self.headers = headers
        self._chunks = chunks
        self._fail_after = fail_after

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        for i, chunk in enumerate(self._chunks):
            if self._fail_after is not None and i >= self._fail_after:
                raise OSError("connection reset")
            yield chunk


def _use(monkeypatch, response):
    monkeypatch.setattr(pipeline_05c.requests, "get",
                        lambda *a, **k: response)


def test_successful_download_is_verified_and_marked(tmp_path, monkeypatch):
    dest = tmp_path / "GSE1" / "GSE1_counts.csv.gz"
    _use(monkeypatch, _Response({"Content-Length": "9"}, [b"abc", b"defghi"]))

    assert pipeline_05c.download_file("ftp://x/GSE1_counts.csv.gz", dest) is True
    assert dest.read_bytes() == b"abcdefghi"
    assert not dest.with_suffix(dest.suffix + ".part").exists()
    assert dest.with_suffix(dest.suffix + ".ok").exists()
    assert pipeline_05c.already_downloaded("GSE1", tmp_path)


def test_mid_stream_failure_leaves_no_file(tmp_path, monkeypatch):
    dest = tmp_path / "GSE1" / "GSE1_counts.csv.gz"
    _use(monkeypatch, _Response({}, [b"abc", b"def"], fail_after=1))

    assert pipeline_05c.download_file("ftp://x/GSE1_counts.csv.gz", dest) is False
    assert not dest.exists()
    assert not dest.with_suffix(dest.suffix + ".part").exists()
    assert not pipeline_05c.already_downloaded("GSE1", tmp_path)


def test_short_read_against_content_length_is_rejected(tmp_path, monkeypatch):
    dest = tmp_path / "GSE1" / "GSE1_counts.csv.gz"
    _use(monkeypatch, _Response({"Content-Length": "100"}, [b"abc"]))

    assert pipeline_05c.download_file("ftp://x/GSE1_counts.csv.gz", dest) is False
    assert not dest.exists()


def test_a_part_file_is_not_counted_as_downloaded(tmp_path):
    suppl = tmp_path / "GSE1"
    suppl.mkdir()
    (suppl / "GSE1_counts.csv.gz.part").write_bytes(b"partial")

    assert not pipeline_05c.already_downloaded("GSE1", tmp_path)


def test_an_unverified_count_file_is_not_counted(tmp_path):
    suppl = tmp_path / "GSE1"
    suppl.mkdir()
    (suppl / "GSE1_counts.csv.gz").write_bytes(b"complete looking")

    assert not pipeline_05c.already_downloaded("GSE1", tmp_path)


def test_a_file_shorter_than_its_marker_is_not_counted(tmp_path):
    suppl = tmp_path / "GSE1"
    suppl.mkdir()
    data = suppl / "GSE1_counts.csv.gz"
    data.write_bytes(b"short")
    pipeline_05c._write_ok_marker(data, 100)

    assert not pipeline_05c.already_downloaded("GSE1", tmp_path)
