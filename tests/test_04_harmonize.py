"""Tests for pipeline/04_harmonize.py and the harmonize package."""

import json

import pytest

from cp_multiomics.harmonize import HARMONIZER_REGISTRY, HarmonizerFactory
from cp_multiomics.harmonize.base import (
    OmicsHarmonizer,
)
from cp_multiomics.harmonize.genomics import GenomicsHarmonizer
from cp_multiomics.harmonize.metabolomics import MetabolomicsHarmonizer
from cp_multiomics.harmonize.transcriptomics import TranscriptomicsHarmonizer

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

def test_all_modalities_registered():
    expected = {"genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"}
    assert expected.issubset(set(HARMONIZER_REGISTRY))


def test_factory_raises_on_unknown():
    with pytest.raises(ValueError, match="Unknown modality"):
        HarmonizerFactory("imaginary")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def test_normalize_platform_rnaseq():
    assert OmicsHarmonizer.normalize_platform("RNA-seq") == "RNA-seq"
    assert OmicsHarmonizer.normalize_platform(
        "Expression profiling by high throughput sequencing") == "RNA-seq"
    assert OmicsHarmonizer.normalize_platform("rnaseq") == "RNA-seq"


def test_normalize_platform_unknown_passthrough():
    result = OmicsHarmonizer.normalize_platform("some novel platform")
    assert result == "some novel platform"


def test_normalize_species_canonical():
    result = OmicsHarmonizer.normalize_species(["homo sapiens", "Mus musculus", "rat"])
    assert "Homo sapiens" in result
    assert "Mus musculus" in result
    assert "Rattus norvegicus" in result


def test_normalize_species_deduplicates():
    result = OmicsHarmonizer.normalize_species(["Homo sapiens", "human", "Homo sapiens"])
    assert result.count("Homo sapiens") == 1


def test_flag_quality_low_n():
    raw = {"n_samples": 2, "year": "2022", "description": "test"}
    flags = OmicsHarmonizer.flag_quality(raw, min_samples=5)
    assert any("low_n_samples" in f for f in flags)


def test_flag_quality_missing_year():
    raw = {"n_samples": 20, "year": "", "description": "test"}
    flags = OmicsHarmonizer.flag_quality(raw, min_samples=5)
    assert "missing_year" in flags


def test_flag_quality_clean_record():
    raw = {"n_samples": 50, "year": "2021", "description": "Chronic pain study"}
    flags = OmicsHarmonizer.flag_quality(raw, min_samples=5)
    assert flags == []


# ---------------------------------------------------------------------------
# Per-modality harmonizers
# ---------------------------------------------------------------------------

CFG = {"filters": {"min_samples": 5}}


def _raw(modality="transcriptomics", **overrides):
    base = {
        "modality": modality,
        "source": "GEO",
        "accession": "GSE99999",
        "title": "Test",
        "description": "Chronic pain",
        "species": ["Homo sapiens"],
        "n_samples": 30,
        "platform": "RNA-seq",
        "year": "2022",
        "url": "https://example.com",
        "keywords": [],
    }
    base.update(overrides)
    return base


def test_transcriptomics_human_only():
    h = TranscriptomicsHarmonizer()
    rec = h._harmonize_one(_raw(), CFG)
    assert rec.has_human
    assert not rec.has_animal
    assert not rec.ortholog_needed
    assert rec.id_space == "HGNC"
    assert rec.effect_size_unit == "log2FC"


def test_transcriptomics_animal_only_needs_ortholog():
    h = TranscriptomicsHarmonizer()
    rec = h._harmonize_one(_raw(species=["Mus musculus"]), CFG)
    assert not rec.has_human
    assert rec.has_animal
    assert rec.ortholog_needed


def test_transcriptomics_both_species_no_ortholog_needed():
    h = TranscriptomicsHarmonizer()
    rec = h._harmonize_one(_raw(species=["Homo sapiens", "Mus musculus"]), CFG)
    assert rec.has_human and rec.has_animal
    # ortholog_needed = has_animal AND NOT has_human → False here
    assert not rec.ortholog_needed


def test_genomics_always_human():
    h = GenomicsHarmonizer()
    rec = h._harmonize_one(
        _raw(modality="genomics", source="GWAS_CATALOG", platform="GWAS", species=[]),
        CFG,
    )
    assert rec.has_human
    assert not rec.has_animal
    assert rec.id_space == "rsID"
    assert not rec.ortholog_needed


def test_genomics_effect_size_or():
    h = GenomicsHarmonizer()
    rec = h._harmonize_one(
        _raw(modality="genomics", source="GWAS_CATALOG",
             description="odds ratio for chronic pain", keywords=[]),
        CFG,
    )
    assert rec.effect_size_unit == "OR"


def test_genomics_small_gwas_flagged():
    h = GenomicsHarmonizer()
    rec = h._harmonize_one(
        _raw(modality="genomics", source="GWAS_CATALOG", n_samples=500),
        CFG,
    )
    assert any("small_gwas" in f for f in rec.quality_flags)


def test_metabolomics_lc_ms_chebi():
    h = MetabolomicsHarmonizer()
    rec = h._harmonize_one(
        _raw(modality="metabolomics", source="METABOLIGHTS", platform="LC-MS"),
        CFG,
    )
    assert rec.id_space == "ChEBI"
    assert rec.effect_size_unit == "SMD"


def test_harmonized_record_jsonl_roundtrip():
    h = TranscriptomicsHarmonizer()
    rec = h._harmonize_one(_raw(), CFG)
    line = rec.to_jsonl()
    d = json.loads(line)
    assert d["accession"] == "GSE99999"
    assert d["id_space"] == "HGNC"
    assert isinstance(d["quality_flags"], list)


# ---------------------------------------------------------------------------
# Pipeline-level (04_harmonize.py)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def pipeline_04(pipeline_03):
    from tests.conftest import load_pipeline_module
    return load_pipeline_module("04_harmonize.py")


def test_load_raw_records_missing_returns_empty(pipeline_04, tmp_path):
    records = pipeline_04.load_raw_records(tmp_path, "transcriptomics")
    assert records == []


def test_load_raw_records_reads_jsonl(pipeline_04, tmp_path):
    mod_dir = tmp_path / "transcriptomics"
    mod_dir.mkdir()
    (mod_dir / "datasets.jsonl").write_text(
        json.dumps({"accession": "GSE1"}) + "\n" +
        json.dumps({"accession": "GSE2"}) + "\n"
    )
    records = pipeline_04.load_raw_records(tmp_path, "transcriptomics")
    assert len(records) == 2


def test_dry_run_returns_zero(pipeline_04, tmp_path):
    cfg = {"filters": {"min_samples": 5}}
    mod_dir = tmp_path / "transcriptomics"
    mod_dir.mkdir()
    (mod_dir / "datasets.jsonl").write_text(
        json.dumps(_raw()) + "\n"
    )
    count = pipeline_04.run_modality(
        "transcriptomics", cfg, tmp_path, tmp_path / "interim", dry_run=True
    )
    assert count == 0
