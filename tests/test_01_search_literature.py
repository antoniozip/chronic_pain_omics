"""Tests for pipeline/01_search_literature.py (logic only, no network calls)."""

import pytest


@pytest.fixture(scope="module")
def m(pipeline_01):
    return pipeline_01


def _make_record(m, source: str, uid: str, doi: str = ""):
    return m.SearchRecord(
        source=source, uid=uid, doi=doi, title="T", authors="A",
        journal="J", year="2024", abstract="", query_label="test",
        retrieved_at="2024-01-01T00:00:00+00:00",
    )


# --- Query building ---

def test_pubmed_query_contains_chronic_pain(m, search_cfg):
    q = m.build_pubmed_query(search_cfg)
    assert "chronic pain" in q


def test_pubmed_query_contains_omics_term(m, search_cfg):
    q = m.build_pubmed_query(search_cfg)
    assert any(term in q for term in ["transcriptomics", "RNA-seq", "proteomics", "metabolomics"])


def test_pubmed_query_contains_date_filter(m, search_cfg):
    q = m.build_pubmed_query(search_cfg)
    assert "PDAT" in q  # date filter injected


def _min_cfg():
    return {
        "query": {
            "terms": ["chronic pain"],
            "omics": {"genomics": "GWAS"},
            "species": {"human": "human"},
        },
        "filters": {"date_range": {"start": "2000-01-01", "end": None}},
    }


def test_topic_query_has_no_date_clause(m):
    q = m.build_topic_query(_min_cfg())
    assert "PDAT" not in q
    assert "FIRST_PDATE" not in q
    assert "chronic pain" in q


def test_pubmed_query_uses_pdat(m):
    q = m.build_query(_min_cfg(), "pubmed")
    assert "[PDAT]" in q


def test_europe_pmc_query_never_uses_pdat(m):
    """[PDAT] silently collapses any Europe PMC query to hitCount=0.

    NOTE: this checks for the bracket-qualified PubMed field tag "[PDAT]"
    rather than the bare substring "PDAT", because "PDAT" is itself a
    substring of "FIRST_PDATE" (F-I-R-S-T-_-P-D-A-T-E) — a bare-substring
    check can never pass once the required FIRST_PDATE clause is present.
    This mirrors the bracket-qualified check already used in
    test_pubmed_query_uses_pdat.
    """
    q = m.build_query(_min_cfg(), "europe_pmc")
    assert "[PDAT]" not in q
    assert "FIRST_PDATE:[2000-01-01 TO " in q


def test_build_pubmed_query_alias_is_preserved(m):
    assert m.build_pubmed_query(_min_cfg()) == m.build_query(_min_cfg(), "pubmed")


def test_unknown_source_raises(m):
    with pytest.raises(ValueError, match="unknown source"):
        m.build_query(_min_cfg(), "scopus")


def test_no_date_filter_when_start_absent(m):
    cfg = _min_cfg()
    cfg["filters"] = {}
    for source in ("pubmed", "europe_pmc"):
        q = m.build_query(cfg, source)
        assert "PDAT" not in q and "FIRST_PDATE" not in q


# --- Deduplication ---

def test_dedup_removes_doi_duplicates(m):
    r1 = _make_record(m, "pubmed", "111", doi="10.1000/xyz")
    r2 = _make_record(m, "europe_pmc", "EPM1", doi="10.1000/xyz")
    unique, n = m.deduplicate([r1, r2])
    assert n == 1 and len(unique) == 1


def test_dedup_keeps_distinct_dois(m):
    r1 = _make_record(m, "pubmed", "111", doi="10.1000/aaa")
    r2 = _make_record(m, "pubmed", "222", doi="10.1000/bbb")
    unique, n = m.deduplicate([r1, r2])
    assert n == 0 and len(unique) == 2


def test_dedup_fallback_uid_when_no_doi(m):
    r1 = _make_record(m, "pubmed", "111")
    r2 = _make_record(m, "pubmed", "111")
    unique, n = m.deduplicate([r1, r2])
    assert n == 1 and len(unique) == 1


def test_dedup_preserves_insertion_order(m):
    records = [_make_record(m, "pubmed", str(i), doi=f"10.1/{i}") for i in range(5)]
    unique, _ = m.deduplicate(records)
    assert [r.uid for r in unique] == [str(i) for i in range(5)]
