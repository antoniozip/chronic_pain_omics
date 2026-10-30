"""Integrity tests for the single-cell stages 16-20.

Covers the review findings that a remote filename could steer a local write
(R-02), a study excluded in step 18 kept a stale groups file (R-03), and an
`--accession` run replaced a whole-corpus artefact with a one-study subset or
truncated it to zero bytes (R-04).
"""

import gzip
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "pipeline" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


step16 = _load("pipeline_16_integrity", "16_download_single_cell.py")
step17 = _load("pipeline_17_integrity", "17_unpack_single_cell.py")
step18 = _load("pipeline_18_integrity", "18_assign_single_cell_groups.py")
step20 = _load("pipeline_20_integrity", "20_harmonize_single_cell.py")


# ---------------------------------------------------------------------------
# R-02 -- remote filenames must not escape the per-series directory
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "../../etc/evil.txt",
    "../evil.txt",
    "/etc/passwd",
    "sub/dir.txt",
    "..",
    ".",
    "",
])
def test_safe_name_rejects_escaping_names(name):
    assert step16.safe_name(name) is None


def test_safe_name_keeps_a_plain_filename():
    assert step16.safe_name("GSE1_RAW.tar") == "GSE1_RAW.tar"


def test_fetch_refuses_a_name_that_escapes_the_target(tmp_path):
    target = tmp_path / "GSE1"
    target.mkdir()
    with pytest.raises(ValueError):
        step16.fetch("https://example/evil", target, "../evil.txt", None)
    assert not (tmp_path / "evil.txt").exists()


# ---------------------------------------------------------------------------
# R-04 -- --accession runs merge into the corpus manifest
# ---------------------------------------------------------------------------

def test_accession_run_merges_the_download_manifest(tmp_path, monkeypatch):
    registry = tmp_path / "registry.csv"
    pd.DataFrame({"accession": ["GSE1", "GSE2"]}).to_csv(registry, index=False)
    out_root = tmp_path / "single_cell"
    out_root.mkdir()
    manifest = out_root / "download_manifest.csv"
    pd.DataFrame([{"accession": "GSE1", "file": "old.tar", "bytes": 1,
                   "outcome": "skipped"}]).to_csv(manifest, index=False)

    monkeypatch.setattr(step16, "REGISTRY", registry)
    monkeypatch.setattr(step16, "OUT_ROOT", out_root)
    monkeypatch.setattr(step16, "MANIFEST", manifest)
    monkeypatch.setattr(step16, "list_files", lambda accession: ["a.txt"])
    monkeypatch.setattr(step16, "remote_size", lambda url: 4)

    def fake_fetch(url, target, name, expected):
        target.mkdir(parents=True, exist_ok=True)
        (target / name).write_bytes(b"1234")
        return "downloaded", 4

    monkeypatch.setattr(step16, "fetch", fake_fetch)
    monkeypatch.setattr("sys.argv", ["16", "--accession", "GSE2"])

    assert step16.main() == 0
    merged = pd.read_csv(manifest)
    assert set(merged["accession"]) == {"GSE1", "GSE2"}


def _write_tar(path: Path, member: str):
    with tarfile.open(path, "w") as tar:
        info = tarfile.TarInfo(member)
        payload = b"counts"
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))


def test_accession_run_merges_the_unpack_manifest(tmp_path, monkeypatch):
    raw_root = tmp_path / "single_cell"
    for acc in ("GSE1", "GSE2"):
        (raw_root / acc).mkdir(parents=True)
        _write_tar(raw_root / acc / f"{acc}_RAW.tar", "GSM9_counts.txt.gz")
    manifest = raw_root / "unpack_manifest.csv"
    pd.DataFrame([{"accession": "GSE1", "gsm": "GSM1", "n_files": 1,
                   "outcome": "unpacked"}]).to_csv(manifest, index=False)

    monkeypatch.setattr(step17, "RAW_ROOT", raw_root)
    monkeypatch.setattr(step17, "MANIFEST", manifest)
    monkeypatch.setattr("sys.argv", ["17", "--accession", "GSE2"])

    assert step17.main() == 0
    merged = pd.read_csv(manifest)
    assert set(merged["accession"]) == {"GSE1", "GSE2"}


# ---------------------------------------------------------------------------
# R-03 -- an excluded study's stale groups file is removed
# ---------------------------------------------------------------------------

def test_excluded_study_loses_its_stale_groups_file(tmp_path):
    groups = tmp_path / "groups"
    groups.mkdir()
    stale = groups / "GSE1_groups.csv"
    stale.write_text("sample_id,group,subject,title\nGSM1,case,GSM1,x\n")

    step18.write_groups("GSE1", [{"sample_id": "GSM1"}], "exclude", groups_dir=groups)
    assert not stale.exists()


def test_included_study_still_gets_a_groups_file(tmp_path):
    groups = tmp_path / "groups"
    groups.mkdir()
    placed = [{"sample_id": "GSM1", "group": "case", "subject": "GSM1", "title": "x"}]

    step18.write_groups("GSE1", placed, "include", groups_dir=groups)

    frame = pd.read_csv(groups / "GSE1_groups.csv")
    assert list(frame["sample_id"]) == ["GSM1"]


def test_accession_run_merges_the_candidate_sheet(tmp_path):
    sheet = tmp_path / "single_cell_candidates.csv"
    pd.DataFrame([{"accession": "GSE1", "verdict": "include"}]).to_csv(sheet, index=False)
    new = pd.DataFrame([{"accession": "GSE2", "verdict": "include"}])

    merged = step18.merge_sheet(sheet, new, replaced={"GSE2"})

    assert set(merged["accession"]) == {"GSE1", "GSE2"}


# ---------------------------------------------------------------------------
# R-04 -- step 20 merges and never writes a zero-byte manifest
# ---------------------------------------------------------------------------

def _harmonize_fixture(tmp_path, monkeypatch, existing):
    registry = tmp_path / "registry.csv"
    pd.DataFrame([
        {"accession": "GSE1", "title": "one", "species": "Homo sapiens",
         "year": 2020, "assay": "10x"},
        {"accession": "GSE2", "title": "two", "species": "Homo sapiens",
         "year": 2021, "assay": "10x"},
    ]).to_csv(registry, index=False)
    candidates = tmp_path / "candidates.csv"
    pd.DataFrame([
        {"accession": "GSE1", "verdict": "include"},
        {"accession": "GSE2", "verdict": "include"},
    ]).to_csv(candidates, index=False)
    interim = tmp_path / "interim"
    pseudobulk = interim / "pseudobulk"
    pseudobulk.mkdir(parents=True)
    manifest = interim / "harmonized.jsonl"

    if existing is not None:
        manifest.write_text("\n".join(json.dumps(r) for r in existing) + "\n")

    monkeypatch.setattr(step20, "REGISTRY", registry)
    monkeypatch.setattr(step20, "CANDIDATES", candidates)
    monkeypatch.setattr(step20, "INTERIM", interim)
    monkeypatch.setattr(step20, "PSEUDOBULK", pseudobulk)
    monkeypatch.setattr(step20, "MANIFEST", manifest)
    return pseudobulk, manifest


def test_accession_run_merges_the_harmonized_manifest(tmp_path, monkeypatch):
    pseudobulk, manifest = _harmonize_fixture(
        tmp_path, monkeypatch, existing=[{"accession": "GSE1", "modality": "single_cell"}])
    with gzip.open(pseudobulk / "GSE2_counts.tsv.gz", "wt") as handle:
        handle.write("feature_id\tS1\nENSG00000000001\t5\n")

    monkeypatch.setattr("sys.argv", ["20", "--accession", "GSE2"])

    assert step20.main() == 0
    records = [json.loads(line) for line in manifest.read_text().splitlines() if line]
    assert {r["accession"] for r in records} == {"GSE1", "GSE2"}


def test_no_records_is_an_error_and_never_a_zero_byte_manifest(tmp_path, monkeypatch):
    _pseudobulk, manifest = _harmonize_fixture(tmp_path, monkeypatch, existing=None)
    monkeypatch.setattr("sys.argv", ["20"])

    assert step20.main() != 0
    assert not manifest.exists()
