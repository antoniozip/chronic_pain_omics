"""Tests for pagination and zero-recall behaviour of the ingester base class."""

from __future__ import annotations

import pytest

from cp_multiomics.ingest.base import DatasetRecord, OmicsIngester


def _record(i: int) -> DatasetRecord:
    return DatasetRecord(
        modality="test", source="TEST", accession=f"ACC{i}", title=f"t{i}",
        description="", species=[], n_samples=1, platform="p", year="2020", url="",
    )


class FakeIngester(OmicsIngester):
    """Serves `total` records in pages, recording how many pages were requested."""

    modality = "test"

    def __init__(self, total: int):
        self.total = total
        self.pages_fetched = 0

    def _build_query(self, terms, species_terms):
        return " OR ".join(terms)

    def _fetch_page(self, query, offset, page_size):
        self.pages_fetched += 1
        remaining = max(0, self.total - offset)
        return [_record(offset + i) for i in range(min(page_size, remaining))]


def test_paginates_past_the_old_500_cap():
    ing = FakeIngester(total=1200)
    out = ing.search_all(["pain"], [], {}, page_size=100)
    assert len(out) == 1200


def test_stops_on_a_short_page():
    ing = FakeIngester(total=250)
    out = ing.search_all(["pain"], [], {}, page_size=100)
    assert len(out) == 250
    assert ing.pages_fetched == 3


def test_hard_limit_truncates_and_warns(caplog):
    ing = FakeIngester(total=1000)
    with caplog.at_level("WARNING"):
        out = ing.search_all(["pain"], [], {}, hard_limit=200, page_size=100)
    assert len(out) == 200
    assert "hard_limit reached" in caplog.text


def test_no_warning_when_hard_limit_not_reached(caplog):
    ing = FakeIngester(total=50)
    with caplog.at_level("WARNING"):
        out = ing.search_all(["pain"], [], {}, hard_limit=200, page_size=100)
    assert len(out) == 50
    assert "hard_limit reached" not in caplog.text


def test_zero_results_raises_by_default():
    ing = FakeIngester(total=0)
    with pytest.raises(RuntimeError, match="returned 0 records"):
        ing.search_all(["pain"], [], {}, page_size=100)


def test_zero_results_allowed_when_opted_in():
    ing = FakeIngester(total=0)
    assert ing.search_all(["pain"], [], {}, page_size=100, allow_empty=True) == []


# --- GenomicsIngester overrides search_all and must honour the same contract ---

def test_genomics_search_all_accepts_hard_limit_keyword():
    """pipeline/03 calls search_all(hard_limit=...); the override must not TypeError."""
    import inspect

    from cp_multiomics.ingest.genomics import GenomicsIngester

    params = inspect.signature(GenomicsIngester.search_all).parameters
    assert "hard_limit" in params
    assert "allow_empty" in params
    assert "max_results" not in params


def test_genomics_search_all_raises_on_zero_recall(monkeypatch):
    from cp_multiomics.ingest.genomics import GenomicsIngester

    ing = GenomicsIngester()
    monkeypatch.setattr(GenomicsIngester, "_get", staticmethod(
        lambda url, params=None, timeout=30: {"_embedded": {"studies": []}}))
    with pytest.raises(RuntimeError, match="returned 0 records"):
        ing.search_all(["chronic pain"], [], {}, page_size=100)


def test_genomics_search_all_warns_when_hard_limit_reached(monkeypatch, caplog):
    from cp_multiomics.ingest.genomics import GenomicsIngester

    ing = GenomicsIngester()
    study = {"accessionId": "GCST1", "diseaseTrait": {"trait": "chronic pain"},
             "publicationInfo": {"publicationDate": "2020-01-01"}}
    monkeypatch.setattr(GenomicsIngester, "_get", staticmethod(
        lambda url, params=None, timeout=30: {"_embedded": {"studies": [study] * 100}}))

    with caplog.at_level("WARNING"):
        out = ing.search_all(["chronic pain"], [], {}, hard_limit=1, page_size=100)
    assert len(out) == 1
    assert "hard_limit reached" in caplog.text


def test_genomics_search_all_warns_when_hard_limit_reached_on_last_term(monkeypatch, caplog):
    """Regression test: the hard_limit warning must fire based on FINAL state, not
    loop position. With PAIN_EFO_TERMS patched to a single term, that term is both
    the first AND the last iteration of the outer `for term in search_terms` loop,
    so the old "warn in the per-term guard before the next iteration" logic had no
    "next" iteration to warn from and silently swallowed the warning."""
    import cp_multiomics.ingest.genomics as genomics_module
    from cp_multiomics.ingest.genomics import GenomicsIngester

    monkeypatch.setattr(genomics_module, "PAIN_EFO_TERMS", ["chronic pain"])

    ing = GenomicsIngester()
    study = {"accessionId": "GCST1", "diseaseTrait": {"trait": "chronic pain"},
             "publicationInfo": {"publicationDate": "2020-01-01"}}
    monkeypatch.setattr(GenomicsIngester, "_get", staticmethod(
        lambda url, params=None, timeout=30: {"_embedded": {"studies": [study] * 100}}))

    with caplog.at_level("WARNING"):
        out = ing.search_all(["chronic pain"], [], {}, hard_limit=1, page_size=100)
    assert len(out) == 1
    assert "hard_limit reached" in caplog.text


def test_genomics_search_all_no_false_positive_warning_on_exhaustion(monkeypatch, caplog):
    """When the source is genuinely exhausted below the cap, no hard_limit warning
    should ever fire — even with a single search term (last-term == only-term)."""
    import cp_multiomics.ingest.genomics as genomics_module
    from cp_multiomics.ingest.genomics import GenomicsIngester

    monkeypatch.setattr(genomics_module, "PAIN_EFO_TERMS", ["chronic pain"])

    ing = GenomicsIngester()
    study = {"accessionId": "GCST1", "diseaseTrait": {"trait": "chronic pain"},
             "publicationInfo": {"publicationDate": "2020-01-01"}}
    # Fewer studies than page_size => short page => genuinely exhausted, well below hard_limit.
    monkeypatch.setattr(GenomicsIngester, "_get", staticmethod(
        lambda url, params=None, timeout=30: {"_embedded": {"studies": [study] * 5}}))

    with caplog.at_level("WARNING"):
        out = ing.search_all(["chronic pain"], [], {}, hard_limit=1000, page_size=100)
    assert len(out) == 1  # deduped by accession
    assert "hard_limit reached" not in caplog.text
