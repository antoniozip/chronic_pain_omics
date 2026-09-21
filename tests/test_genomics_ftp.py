import pytest

from cp_multiomics.genomics.ftp import (
    FileTooLargeError,
    download_file,
    find_harmonised_file,
    harmonised_dir,
)

FTP = "https://ftp.example/summary_statistics"


def test_harmonised_dir_range_math():
    # GCST90091913 -> range GCST90091001-GCST90092000
    d = harmonised_dir("GCST90091913", FTP)
    assert d == f"{FTP}/GCST90091001-GCST90092000/GCST90091913/harmonised/"


def test_harmonised_dir_low_accession():
    # GCST006927 -> GCST006001-GCST007000
    d = harmonised_dir("GCST006927", FTP)
    assert d == f"{FTP}/GCST006001-GCST007000/GCST006927/harmonised/"


def test_find_prefers_canonical_name():
    files = ["GCST90091913.h.tsv.gz", "34737426-GCST90091913-EFO_1-Build37.h.tsv.gz"]
    url = find_harmonised_file("GCST90091913", FTP, lambda u: files)
    assert url.endswith("/GCST90091913.h.tsv.gz")


def test_find_falls_back_to_any_harmonised():
    files = ["34737426-GCST90043740-HP_0012532.h.tsv.gz", "README.txt"]
    url = find_harmonised_file("GCST90043740", FTP, lambda u: files)
    assert url.endswith("34737426-GCST90043740-HP_0012532.h.tsv.gz")


def test_find_returns_none_when_no_harmonised_file():
    assert find_harmonised_file("GCST006927", FTP, lambda u: []) is None
    assert find_harmonised_file("GCST006927", FTP, lambda u: ["README.txt"]) is None


def test_download_writes_bytes(tmp_path):
    dest = tmp_path / "x.tsv.gz"
    def stream(url):
        return {"Content-Length": "9"}, [b"abc", b"def", b"ghi"]
    n = download_file("http://x", dest, stream, max_bytes=1_000_000)
    assert n == 9
    assert dest.read_bytes() == b"abcdefghi"


def test_download_aborts_when_too_large(tmp_path):
    dest = tmp_path / "big.tsv.gz"
    def stream(url):
        return {"Content-Length": str(2_000_000_000)}, [b"x"]
    with pytest.raises(FileTooLargeError):
        download_file("http://x", dest, stream, max_bytes=1_500_000_000)
    assert not dest.exists()  # nothing written on abort
