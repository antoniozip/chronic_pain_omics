"""Tests for pipeline/03_ingest_omics.py and the ingest package (no network calls)."""

import json

import pytest

from cp_multiomics.ingest import INGESTER_REGISTRY, IngesterFactory
from cp_multiomics.ingest.base import DatasetRecord
from cp_multiomics.ingest.genomics import GenomicsIngester
from cp_multiomics.ingest.metabolomics import MetabolomicsIngester
from cp_multiomics.ingest.transcriptomics import TranscriptomicsIngester

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

def test_all_modalities_registered():
    expected = {"genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"}
    assert expected.issubset(set(INGESTER_REGISTRY))


def test_factory_returns_correct_type():
    ingester = IngesterFactory("transcriptomics")
    assert isinstance(ingester, TranscriptomicsIngester)


def test_factory_raises_on_unknown_modality():
    with pytest.raises(ValueError, match="Unknown modality"):
        IngesterFactory("unknown_modality")


# ---------------------------------------------------------------------------
# DatasetRecord
# ---------------------------------------------------------------------------

def _make_record(**kwargs) -> DatasetRecord:
    defaults = dict(
        modality="transcriptomics", source="GEO", accession="GSE99999",
        title="Test study", description="", species=["Homo sapiens"],
        n_samples=20, platform="RNA-seq", year="2022", url="https://example.com",
    )
    defaults.update(kwargs)
    return DatasetRecord(**defaults)


def test_dataset_record_to_jsonl_roundtrip():
    r = _make_record()
    line = r.to_jsonl()
    d = json.loads(line)
    assert d["accession"] == "GSE99999"
    assert d["n_samples"] == 20
    assert d["species"] == ["Homo sapiens"]


def test_dataset_record_has_retrieved_at():
    r = _make_record()
    assert r.retrieved_at  # auto-populated


# ---------------------------------------------------------------------------
# Query builders
# ---------------------------------------------------------------------------

def test_transcriptomics_query_contains_pain_and_rna():
    ingester = TranscriptomicsIngester()
    q = ingester._build_query(["chronic pain"], ["Homo sapiens"])
    assert "chronic pain" in q
    assert "Homo sapiens" in q
    assert "Expression profiling" in q or "gse" in q


def test_genomics_query_joins_terms():
    ingester = GenomicsIngester()
    q = ingester._build_query(["chronic pain", "neuropathic pain"], ["Homo sapiens"])
    assert "chronic pain" in q


def test_metabolomics_excludes_lipidomics():
    ingester = MetabolomicsIngester()
    assert ingester._is_lipidomics("this study is about lipidomics profiling")
    assert not ingester._is_lipidomics("metabolomics in chronic pain patients")


# ---------------------------------------------------------------------------
# Filters (pipeline level)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def pipeline_03(pipeline_02):
    from tests.conftest import load_pipeline_module
    return load_pipeline_module("03_ingest_omics.py")


def test_apply_filters_drops_small_studies(pipeline_03):
    records = [
        _make_record(n_samples=3, year="2020"),
        _make_record(n_samples=50, year="2020"),
    ]
    kept, dropped = pipeline_03.apply_filters(records, {"min_samples": 5, "year_min": 2000})
    assert dropped == 1
    assert len(kept) == 1


def test_apply_filters_drops_old_studies(pipeline_03):
    records = [
        _make_record(n_samples=10, year="1998"),
        _make_record(n_samples=10, year="2010"),
    ]
    kept, dropped = pipeline_03.apply_filters(records, {"min_samples": 5, "year_min": 2000})
    assert dropped == 1 and len(kept) == 1


def test_apply_filters_handles_missing_year(pipeline_03):
    r = _make_record(n_samples=10, year="")
    kept, dropped = pipeline_03.apply_filters([r], {"min_samples": 5, "year_min": 2000})
    # No year → year=0 → not penalized (we don't know; include by default)
    assert len(kept) == 1


# ---------------------------------------------------------------------------
# Manifest I/O
# ---------------------------------------------------------------------------

def test_write_manifest_creates_jsonl(pipeline_03, tmp_path):
    records = [_make_record(accession=f"GSE{i}") for i in range(3)]
    out = tmp_path / "transcriptomics" / "datasets.jsonl"
    pipeline_03.write_manifest(records, out)
    assert out.exists()
    lines = out.read_text().strip().split("\n")
    assert len(lines) == 3
    assert json.loads(lines[0])["accession"] == "GSE0"


# ---------------------------------------------------------------------------
# GWAS sample count parser
# ---------------------------------------------------------------------------

def test_gwas_parse_sample_count():
    ingester = GenomicsIngester()
    assert ingester._parse_sample_count("4,326 European ancestry") == 4326
    assert ingester._parse_sample_count("") == 0
    assert ingester._parse_sample_count("N/A") == 0


# ---------------------------------------------------------------------------
# main() resilience: one modality's zero-recall RuntimeError must not abort
# the rest, but must still surface as a non-zero exit code.
# ---------------------------------------------------------------------------

def test_main_continues_past_one_failed_modality_and_exits_nonzero(
    pipeline_03, monkeypatch, tmp_path
):
    calls = []

    def fake_run_modality(modality, cfg, dry_run):
        calls.append(modality)
        if modality == "metabolomics":
            raise RuntimeError(f"[{modality}] search returned 0 records for terms [...]")
        return 1

    monkeypatch.setattr(pipeline_03, "run_modality", fake_run_modality)
    monkeypatch.setattr(
        pipeline_03, "load_config", lambda path: {"search": {}, "filters": {}, "output": {}}
    )
    monkeypatch.setattr(
        "sys.argv", ["03_ingest_omics.py", "--modality", "all"]
    )

    with pytest.raises(SystemExit) as exc_info:
        pipeline_03.main()

    assert exc_info.value.code == 1
    # Every modality was attempted despite metabolomics failing.
    assert set(calls) == set(pipeline_03.ALL_MODALITIES)
