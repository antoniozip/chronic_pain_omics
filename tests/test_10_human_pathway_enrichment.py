"""An Enrichr outage must not be written out as a three-way comparison.

`build_threeway_comparison` reads "the library returned nothing significant"
as "the arm found nothing", so a failed `concordant` call produced a table in
which every term was `in_concordant=False` — an authoritative-looking record
that the translational arm found nothing, which is the manuscript's central
claim. Error and no-overlap are different facts and must stay distinguishable.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tests.conftest import load_pipeline_module


@pytest.fixture(scope="module")
def pipeline_10():
    return load_pipeline_module("10_human_pathway_enrichment.py")


def _terms(term: str = "T1", padj: float = 0.01) -> pd.DataFrame:
    return pd.DataFrame({
        "library": ["KEGG_2021_Human"],
        "term": [term],
        "adjusted_pval": [padj],
        "combined_score": [10.0],
        "n_genes": [3],
        "overlapping_genes": [["A", "B", "C"]],
    })


def test_builder_refuses_to_write_when_an_arm_is_unavailable(pipeline_10, tmp_path):
    out = tmp_path / "threeway_comparison.csv"

    with pytest.raises(ValueError, match="concordant"):
        pipeline_10.build_threeway_comparison(
            [_terms()], [_terms()], [], out_path=out,
            unavailable={"concordant": "enrichr_error"},
        )

    assert not out.exists(), "a failed arm must not leave a comparison table"


def test_no_terms_from_a_working_arm_is_written_as_no_overlap(pipeline_10, tmp_path):
    """An empty result list is a valid "nothing survived", not a failure."""
    out = tmp_path / "threeway_comparison.csv"

    df = pipeline_10.build_threeway_comparison(
        [_terms()], [], [], out_path=out,
    )

    assert out.exists()
    assert not df["in_concordant"].any()


def test_enrichr_error_aborts_the_arm_instead_of_reporting_no_overlap(
        pipeline_10, tmp_path, monkeypatch):
    def fake(gene_list, library, list_name="query"):
        if library == "FAILED_LIB":
            return None
        return _terms()

    monkeypatch.setattr(pipeline_10, "enrichr_run", fake)
    monkeypatch.setattr(pipeline_10.time, "sleep", lambda _s: None)

    with pytest.raises(SystemExit, match="Refusing"):
        pipeline_10.run_enrichr_arm(
            "concordant", list("ABCDE"), "concordant",
            libs=["FAILED_LIB", "KEGG_2021_Human"], out_dir=tmp_path,
        )

    assert not any(tmp_path.iterdir())


def test_empty_result_from_a_working_arm_does_not_abort(pipeline_10, tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline_10, "enrichr_run",
                        lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(pipeline_10.time, "sleep", lambda _s: None)

    assert pipeline_10.run_enrichr_arm(
        "concordant", list("ABCDE"), "concordant",
        libs=["KEGG_2021_Human"], out_dir=tmp_path,
    ) == []
