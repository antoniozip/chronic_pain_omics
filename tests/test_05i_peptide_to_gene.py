"""Tests for pipeline/05i_peptide_to_gene.py (R-06).

`06_meta_analysis.R` prefers `*_effects_gene.csv` over the other per-study
tables, so a run that maps nothing must not exit 0 while leaving the previous
run's gene table in place to be consumed as if it were fresh.
"""

import importlib.util
from pathlib import Path

import pandas as pd

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_05i",
    Path(__file__).parent.parent / "pipeline" / "05i_peptide_to_gene.py",
)
pipeline_05i = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_05i)


def _paths(tmp_path, monkeypatch, map_rows, input_rows):
    mapping = tmp_path / "peptide_gene_map.csv"
    pd.DataFrame(map_rows, columns=["peptide", "gene_symbol"]).to_csv(mapping, index=False)
    source = tmp_path / "PXD013362_effects_collapsed.csv"
    pd.DataFrame(input_rows, columns=["feature_id", "pval", "study_id"]).to_csv(
        source, index=False)
    out = tmp_path / "PXD013362_effects_gene.csv"
    monkeypatch.setattr(pipeline_05i, "MAP", mapping)
    monkeypatch.setattr(pipeline_05i, "IN", source)
    monkeypatch.setattr(pipeline_05i, "OUT", out)
    return out


def test_no_assignments_removes_the_stale_gene_table(tmp_path, monkeypatch):
    out = _paths(
        tmp_path, monkeypatch,
        map_rows=[("", "CALCA"), ("", "PENK")],
        input_rows=[("PEPTIDE1", "0.01", "unit")],
    )
    out.write_text("feature_id,effect_size,se,pval,study_id\nSTALE,1,1,0.5,unit\n")

    assert pipeline_05i.main() != 0
    assert not out.exists()


def test_no_match_removes_the_stale_gene_table(tmp_path, monkeypatch):
    out = _paths(
        tmp_path, monkeypatch,
        map_rows=[("PEPTIDE1", "CALCA")],
        input_rows=[("OTHER", "0.01", "unit")],
    )
    out.write_text("feature_id,effect_size,se,pval,study_id\nSTALE,1,1,0.5,unit\n")

    assert pipeline_05i.main() != 0
    assert not out.exists()


def test_a_successful_mapping_still_writes_the_gene_table(tmp_path, monkeypatch):
    out = _paths(
        tmp_path, monkeypatch,
        map_rows=[("PEPTIDE1", "CALCA")],
        input_rows=[("PEPTIDE1", "0.01", "unit")],
    )

    assert pipeline_05i.main() == 0
    assert out.exists()
    assert list(pd.read_csv(out)["feature_id"]) == ["CALCA"]
