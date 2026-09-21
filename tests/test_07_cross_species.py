"""Tests for ortholog package and pipeline/07_cross_species.py."""

import json
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from cp_multiomics.ortholog.concordance import build_concordance_table, compute_summary
from cp_multiomics.ortholog.hcop import OrthologMapper, OrthologRecord

# ---------------------------------------------------------------------------
# OrthologMapper — cache I/O
# ---------------------------------------------------------------------------

def test_mapper_loads_empty_cache(tmp_path):
    mapper = OrthologMapper(tmp_path / "cache.json")
    assert mapper._cache == {}


def test_mapper_saves_and_reloads_cache(tmp_path):
    cache_path = tmp_path / "cache.json"
    mapper = OrthologMapper(cache_path)

    # Inject a record directly
    mapper._cache["Ptgs2:10090"] = OrthologRecord(
        source_symbol="Ptgs2", source_species="Mus musculus",
        human_symbol="PTGS2", human_entrez=5743,
        confidence="high", sources=["PTGS2"],
    )
    mapper._save_cache()

    mapper2 = OrthologMapper(cache_path)
    assert "Ptgs2:10090" in mapper2._cache
    assert mapper2._cache["Ptgs2:10090"].human_symbol == "PTGS2"


def test_mapper_raises_on_unsupported_species(tmp_path):
    mapper = OrthologMapper(tmp_path / "cache.json")
    with pytest.raises(ValueError, match="Unsupported species"):
        mapper.map_to_human(["Gene1"], from_species="Drosophila melanogaster")


def test_mapper_uses_cache_without_network(tmp_path):
    cache_path = tmp_path / "cache.json"
    mapper = OrthologMapper(cache_path)
    mapper._cache["Tnf:10090"] = OrthologRecord(
        source_symbol="Tnf", source_species="Mus musculus",
        human_symbol="TNF", human_entrez=7124,
        confidence="high", sources=["TNF"],
    )
    mapper._save_cache()

    # Fresh mapper should serve from cache, no network call
    mapper2 = OrthologMapper(cache_path)
    with patch.object(mapper2, "_fetch_and_cache") as mock_fetch:
        result = mapper2.map_to_human(["Tnf"], from_species="Mus musculus")
        mock_fetch.assert_not_called()

    assert result["Tnf"] is not None
    assert result["Tnf"].human_symbol == "TNF"


def test_mapper_returns_none_for_missing_symbol(tmp_path):
    mapper = OrthologMapper(tmp_path / "cache.json")
    mapper._cache["NonExistent:10090"] = None
    mapper._save_cache()

    mapper2 = OrthologMapper(tmp_path / "cache.json")
    with patch.object(mapper2, "_fetch_and_cache"):
        result = mapper2.map_to_human(["NonExistent"], from_species="Mus musculus")
    assert result["NonExistent"] is None


# ---------------------------------------------------------------------------
# concordance.build_concordance_table
# ---------------------------------------------------------------------------

def _pooled(feature_ids, yi_values, padj_values, source="human"):
    return pd.DataFrame({
        "feature_id": feature_ids,
        "yi_pooled":  yi_values,
        "padj_pooled": padj_values,
        "k": [3] * len(feature_ids),
    })


def _make_mapper_with_map(mapping: dict, tmp_path) -> OrthologMapper:
    """Return OrthologMapper whose map_to_human returns a fixed mapping."""
    mapper = OrthologMapper(tmp_path / "cache.json")
    records = {
        sym: OrthologRecord(sym, "Mus musculus", human, None, "high", [human])
        for sym, human in mapping.items()
    }
    mapper.map_to_human = MagicMock(return_value=records)
    return mapper


def test_concordance_table_same_direction(tmp_path):
    human  = _pooled(["GENE1", "GENE2"], [1.0, -0.5], [0.01, 0.03])
    animal = _pooled(["Gene1", "Gene2"], [0.8, -0.4], [0.02, 0.04])
    mapper = _make_mapper_with_map({"Gene1": "GENE1", "Gene2": "GENE2"}, tmp_path)

    result = build_concordance_table(human, animal, "Mus musculus", mapper, padj_threshold=0.05)
    assert len(result) == 2
    assert result["concordant"].all()


def test_concordance_table_opposite_direction_flagged(tmp_path):
    human  = _pooled(["GENE1"], [1.0],  [0.01])
    animal = _pooled(["Gene1"], [-0.8], [0.02])
    mapper = _make_mapper_with_map({"Gene1": "GENE1"}, tmp_path)

    result = build_concordance_table(human, animal, "Mus musculus", mapper)
    assert len(result) == 1
    assert not result["concordant"].iloc[0]


def test_concordance_table_empty_when_no_overlap(tmp_path):
    human  = _pooled(["GENE_X"], [1.0], [0.01])
    animal = _pooled(["GeneY"],  [0.9], [0.02])
    mapper = _make_mapper_with_map({"GeneY": "GENE_Y_DIFFERENT"}, tmp_path)

    result = build_concordance_table(human, animal, "Mus musculus", mapper)
    assert result.empty


def test_concordance_table_empty_animal_df_returns_empty(tmp_path):
    human  = _pooled(["GENE1"], [1.0], [0.01])
    animal = pd.DataFrame()
    mapper = OrthologMapper(tmp_path / "cache.json")

    result = build_concordance_table(human, animal, "Mus musculus", mapper)
    assert result.empty


# ---------------------------------------------------------------------------
# concordance.compute_summary
# ---------------------------------------------------------------------------

def test_compute_summary_pct_concordant():
    df = pd.DataFrame({
        "human_yi":  [1.0, -0.5,  0.8, -0.3],
        "animal_yi": [0.8, -0.4, -0.6,  0.2],
        "concordant": [True, True, False, False],
        "both_significant": [True, True, False, False],
        "human_padj": [0.01, 0.01, 0.1, 0.2],
        "animal_padj": [0.02, 0.02, 0.3, 0.4],
    })
    summary = compute_summary(df, "transcriptomics")
    assert summary.pct_concordant == pytest.approx(50.0)
    assert summary.n_ortholog_pairs == 4
    assert summary.n_concordant == 2
    assert not any(v != v for v in [summary.pearson_r, summary.spearman_r])  # not NaN


def test_compute_summary_handles_too_few_pairs():
    df = pd.DataFrame({"human_yi": [1.0], "animal_yi": [0.8],
                       "concordant": [True], "both_significant": [True],
                       "human_padj": [0.01], "animal_padj": [0.02]})
    summary = compute_summary(df, "transcriptomics")
    assert summary.n_ortholog_pairs == 1
    import math
    assert math.isnan(summary.pearson_r)


# ---------------------------------------------------------------------------
# split_by_species — species attribution of pooled rows
# ---------------------------------------------------------------------------

def _write_manifest(tmp_path, modality="transcriptomics"):
    """Manifest with one human-only study, one mouse study and one rat study."""
    mod_dir = tmp_path / modality
    mod_dir.mkdir(parents=True, exist_ok=True)
    rows = [
        {"accession": "HUM1", "has_human": True,  "has_animal": False,
         "species_canonical": ["Homo sapiens"]},
        {"accession": "MOU1", "has_human": False, "has_animal": True,
         "species_canonical": ["Mus musculus"]},
        {"accession": "RAT1", "has_human": False, "has_animal": True,
         "species_canonical": ["Rattus norvegicus"]},
    ]
    with open(mod_dir / "harmonized.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    return tmp_path


def _pooled_with_studies(features, study_ids):
    return pd.DataFrame({
        "feature_id": features,
        "study_ids": study_ids,
        "yi_pooled": [1.0] * len(features),
        "padj_pooled": [0.01] * len(features),
        "k": [len(s.split(";")) for s in study_ids],
    })


def test_split_by_species_sides_are_disjoint(tmp_path, pipeline_07):
    """A pooled row must never land on both sides of the comparison.

    A feature pooled across human and mouse studies is a single, species-
    agnostic estimate. Placing it in both subsets compares it with itself,
    which is concordant by construction.
    """
    _write_manifest(tmp_path)
    pooled = _pooled_with_studies(
        ["MIXED", "HUMAN_ONLY", "Mouse_only"],
        ["HUM1;MOU1", "HUM1", "MOU1"],
    )

    human_df, animal_df = pipeline_07.split_by_species(
        pooled, tmp_path, "transcriptomics", "Mus musculus"
    )

    overlap = set(human_df["feature_id"]) & set(animal_df["feature_id"])
    assert overlap == set(), f"feature(s) on both sides: {overlap}"


def test_split_by_species_excludes_cross_species_rows(tmp_path, pipeline_07):
    """A row pooling both species belongs to neither single-species subset."""
    _write_manifest(tmp_path)
    pooled = _pooled_with_studies(
        ["MIXED", "HUMAN_ONLY", "Mouse_only"],
        ["HUM1;MOU1", "HUM1", "MOU1"],
    )

    human_df, animal_df = pipeline_07.split_by_species(
        pooled, tmp_path, "transcriptomics", "Mus musculus"
    )

    assert list(human_df["feature_id"]) == ["HUMAN_ONLY"]
    assert list(animal_df["feature_id"]) == ["Mouse_only"]


def test_split_by_species_human_side_excludes_any_animal_species(tmp_path, pipeline_07):
    """The human side must exclude rows drawing on a non-target animal too.

    Filtering only against the target species would let a human+rat row count
    as human when the target species is mouse.
    """
    _write_manifest(tmp_path)
    pooled = _pooled_with_studies(["HUM_PLUS_RAT"], ["HUM1;RAT1"])

    human_df, _ = pipeline_07.split_by_species(
        pooled, tmp_path, "transcriptomics", "Mus musculus"
    )

    assert human_df.empty


# ---------------------------------------------------------------------------
# load_species_pooled — species-stratified pooled effects (06g)
# ---------------------------------------------------------------------------

def _write_species_pooled(meta_dir, modality, slug, features):
    d = meta_dir / modality / "by_species"
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "feature_id": features,
        "k": [3] * len(features),
        "study_ids": ["A;B;C"] * len(features),
        "yi": [0.5] * len(features),
        "se": [0.1] * len(features),
        "pval": [0.001] * len(features),
        "padj": [0.01] * len(features),
        "species": [slug.replace("_", " ")] * len(features),
    }).to_csv(d / f"{slug}_pooled.csv", index=False)


def test_load_species_pooled_renames_to_pooled_schema(tmp_path, pipeline_07):
    """06g emits yi/padj; the concordance builder consumes yi_pooled/padj_pooled."""
    _write_species_pooled(tmp_path, "transcriptomics", "Mus_musculus", ["Gene1", "Gene2"])

    df = pipeline_07.load_species_pooled(tmp_path, "transcriptomics", "Mus musculus")

    assert df is not None
    assert {"feature_id", "yi_pooled", "padj_pooled", "k"} <= set(df.columns)
    assert len(df) == 2


def test_load_species_pooled_returns_none_when_absent(tmp_path, pipeline_07):
    (tmp_path / "transcriptomics").mkdir(parents=True)
    assert pipeline_07.load_species_pooled(
        tmp_path, "transcriptomics", "Mus musculus"
    ) is None


# ---------------------------------------------------------------------------
# ortholog_mapping_stats — coverage of the ortholog map
# ---------------------------------------------------------------------------

def test_mapping_stats_counts_symbols_with_a_human_ortholog(tmp_path, pipeline_07):
    animal = _pooled(["Gene1", "Gene2", "Probe3"], [1.0, -0.5, 0.2], [0.01, 0.03, 0.5])
    mapper = _make_mapper_with_map({"Gene1": "GENE1", "Gene2": "GENE2"}, tmp_path)

    stats = pipeline_07.ortholog_mapping_stats(animal, "Mus musculus", mapper)

    assert stats["n_animal_features"] == 3
    assert stats["n_mapped_to_human"] == 2
    assert stats["pct_mapped_to_human"] == pytest.approx(200 / 3)


def test_mapping_stats_ignores_repeated_symbols(tmp_path, pipeline_07):
    animal = _pooled(["Gene1", "Gene1", "Probe3"], [1.0, 1.0, 0.2], [0.01, 0.01, 0.5])
    mapper = _make_mapper_with_map({"Gene1": "GENE1"}, tmp_path)

    stats = pipeline_07.ortholog_mapping_stats(animal, "Mus musculus", mapper)

    assert stats["n_animal_features"] == 2
    assert stats["n_mapped_to_human"] == 1


def test_mapping_stats_on_empty_input(tmp_path, pipeline_07):
    mapper = OrthologMapper(tmp_path / "cache.json")
    stats = pipeline_07.ortholog_mapping_stats(pd.DataFrame(), "Mus musculus", mapper)
    assert stats["n_animal_features"] == 0
    assert stats["n_mapped_to_human"] == 0
    import math
    assert math.isnan(stats["pct_mapped_to_human"])
