"""The PRIDE v2 keyword parameter is ignored; discovery must go via EBI Search."""

from __future__ import annotations

import logging

import pytest

from cp_multiomics.ingest.proteomics import ProteomicsIngester


@pytest.fixture
def ingester(monkeypatch):
    ing = ProteomicsIngester()

    search_payload = {
        "hitCount": 2,
        "entries": [
            {"id": "PXD054342", "source": "pride"},
            {"id": "PXD013362", "source": "pride"},
        ],
    }
    metadata = {
        "PXD054342": {
            "accession": "PXD054342",
            "title": "Proteomic Analysis of Spinal Dorsal Horn in Prior-Exercise",
            "projectDescription": "Neuropathic pain study.",
            "organisms": [{"name": "Mus musculus (mouse)"}],
            "instruments": [{"name": "LTQ Orbitrap Elite"}],
            "publicationDate": "2025-05-07",
            "samplesCount": None,
            "keywords": ["neuropathic pain"],
        },
        "PXD013362": {
            "accession": "PXD013362",
            "title": "PACAP and other neuropeptide targets link chronic migraine",
            "projectDescription": "Migraine and OIH.",
            "organisms": [{"name": "Mus musculus (mouse)"}],
            "instruments": [{"name": "LTQ Orbitrap"}],
            "publicationDate": "2019-11-06",
            "samplesCount": 12,
            "keywords": [],
        },
    }

    def fake_get(url, params=None, timeout=30):
        if "ebisearch" in url:
            return search_payload
        acc = url.rstrip("/").split("/")[-1]
        return metadata[acc]

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake_get))
    return ing


def test_discovery_uses_ebi_search_not_the_ignored_pride_keyword(ingester, monkeypatch):
    seen: list[str] = []
    original = ProteomicsIngester._get

    def spy(url, params=None, timeout=30):
        seen.append(url)
        return original(url, params, timeout)

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(spy))
    ingester._fetch_page('"chronic pain"', 0, 10)

    assert any("ebisearch" in u for u in seen)
    assert not any(u.endswith("/archive/v2/projects") for u in seen)


def test_records_carry_real_titles(ingester):
    records = ingester._fetch_page('"chronic pain"', 0, 10)
    assert len(records) == 2
    titles = {r.title for r in records}
    assert "PACAP and other neuropeptide targets link chronic migraine" in titles
    assert all(r.title for r in records), "empty titles mean EBI Search fields were trusted"


def test_null_sample_count_becomes_minus_one(ingester):
    records = {r.accession: r for r in ingester._fetch_page('"chronic pain"', 0, 10)}
    assert records["PXD054342"].n_samples == -1   # samplesCount was None
    assert records["PXD013362"].n_samples == 12


def test_species_and_platform_are_populated(ingester):
    rec = next(r for r in ingester._fetch_page('"chronic pain"', 0, 10)
               if r.accession == "PXD054342")
    assert rec.species == ["Mus musculus (mouse)"]
    assert rec.platform == "LTQ Orbitrap Elite"
    assert rec.source == "PRIDE"
    assert rec.year == "2025"


def test_empty_search_returns_empty_list(ingester, monkeypatch):
    def empty(url, params=None, timeout=30):
        return {"entries": []}

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(empty))
    assert ingester._fetch_page('"chronic pain"', 0, 10) == []


def test_metadata_failure_skips_that_accession_only(ingester, monkeypatch):
    def flaky(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXD054342"}, {"id": "PXDBROKEN"}]}
        if url.endswith("PXDBROKEN"):
            raise RuntimeError("404")
        return {
            "accession": "PXD054342", "title": "Spinal Dorsal Horn",
            "projectDescription": "", "organisms": [], "instruments": [],
            "publicationDate": "2025-05-07", "samplesCount": 3, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(flaky))
    records = ingester._fetch_page('"chronic pain"', 0, 10)
    assert [r.accession for r in records] == ["PXD054342"]


# ---------------------------------------------------------------------------
# Hazard tests: base.py's search_all() treats a page SHORTER than requested
# as "source exhausted" and stops. A single transient metadata failure must
# not be able to silently truncate the whole ingestion run.
# ---------------------------------------------------------------------------


def test_metadata_failure_is_retried_once_before_being_dropped(ingester, monkeypatch):
    """A transient failure that succeeds on retry must not drop the record."""
    calls = {"PXDRETRY": 0}

    def flaky_once(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXDRETRY"}]}
        calls["PXDRETRY"] += 1
        if calls["PXDRETRY"] == 1:
            raise RuntimeError("transient 503")
        return {
            "accession": "PXDRETRY", "title": "Recovered after retry",
            "projectDescription": "", "organisms": [], "instruments": [],
            "publicationDate": "2024-01-01", "samplesCount": 5, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(flaky_once))
    records = ingester._fetch_page('"chronic pain"', 0, 10)

    assert calls["PXDRETRY"] == 2, "must retry exactly once before giving up"
    assert [r.accession for r in records] == ["PXDRETRY"]


def test_persistent_metadata_failure_logs_metadata_dropped_warning(ingester, monkeypatch, caplog):
    """A failure that persists through the retry must be logged, not silently swallowed."""

    def always_fails(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXDBROKEN"}]}
        raise RuntimeError("404")

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(always_fails))
    with caplog.at_level(logging.WARNING):
        records = ingester._fetch_page('"chronic pain"', 0, 10)

    assert records == []
    assert any(
        "metadata dropped" in rec.message and "PXDBROKEN" in rec.message
        for rec in caplog.records
    ), "dropped accessions must be logged with the literal substring 'metadata dropped'"


def test_page_shrink_relative_to_discovered_accessions_is_logged(ingester, monkeypatch, caplog):
    """If the page returns fewer records than accessions discovered, log 'page shrank'
    with both counts so search_all's exhaustion heuristic never truncates silently.
    """

    def one_broken_of_two(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXD054342"}, {"id": "PXDBROKEN"}]}
        if url.endswith("PXDBROKEN"):
            raise RuntimeError("404")
        return {
            "accession": "PXD054342", "title": "Spinal Dorsal Horn",
            "projectDescription": "", "organisms": [], "instruments": [],
            "publicationDate": "2025-05-07", "samplesCount": 3, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(one_broken_of_two))
    with caplog.at_level(logging.WARNING):
        records = ingester._fetch_page('"chronic pain"', 0, 10)

    assert len(records) == 1
    assert any(
        "page shrank" in rec.message and "2" in rec.message and "1" in rec.message
        for rec in caplog.records
    ), "a shrunk page must be logged with 'page shrank' and both the discovered/built counts"


# ---------------------------------------------------------------------------
# Backfill regression tests: one dropped accession must not truncate the
# whole multi-page run (the reviewer's reproduction: 250 accessions, one
# permanently 404-ing, yielded only 99 records instead of 249).
# ---------------------------------------------------------------------------


def _paginated_corpus_get(corpus: list[str], broken: set[str]):
    """Build a fake `_get` that serves EBI Search pages honoring start/size,
    and PRIDE metadata for any accession not in `broken` (which always 404s).
    """

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            start = params["start"]
            size = params["size"]
            chunk = corpus[start:start + size]
            return {"entries": [{"id": a} for a in chunk]}
        acc = url.rstrip("/").split("/")[-1]
        if acc in broken:
            raise RuntimeError("404")
        return {
            "accession": acc, "title": f"Study {acc}", "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2022-01-01",
            "samplesCount": 4, "keywords": [],
        }

    return fake


def test_one_dropped_accession_does_not_truncate_a_250_study_corpus(ingester, monkeypatch):
    """Reviewer's exact reproduction: 250 accessions, page_size=100, the 42nd
    accession permanently fails metadata. Must yield 249 records, not 99.
    """
    corpus = [f"PXD1{i:06d}" for i in range(250)]
    broken = {corpus[41]}  # the 42nd accession

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )
    records = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(records) == 249
    assert corpus[41] not in {r.accession for r in records}


def test_fetch_page_dedups_via_seen_accessions(ingester, monkeypatch):
    """An accession already emitted must never be emitted again, since
    backfill can consume accessions beyond the caller's requested window and
    search_all advances its own offset by len(batch)."""

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXD1"}, {"id": "PXD2"}]}
        acc = url.rstrip("/").split("/")[-1]
        return {
            "accession": acc, "title": acc, "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2020-01-01",
            "samplesCount": 1, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))

    first = ingester._fetch_page('"q"', 0, 2)
    assert {r.accession for r in first} == {"PXD1", "PXD2"}

    second = ingester._fetch_page('"q"', 0, 2)
    assert second == [], "both accessions were already seen; nothing new to emit"


def test_genuine_exhaustion_still_terminates_without_looping(ingester, monkeypatch):
    """A corpus smaller than page_size must return exactly what exists and
    stop, not loop waiting for records that will never arrive."""
    corpus = [f"PXD2{i:06d}" for i in range(30)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, set()))
    )
    records = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(records) == 30


def test_backfill_exhausted_warning_and_termination(ingester, monkeypatch, caplog):
    """If every accession keeps failing metadata, the backfill loop must
    still terminate (bounded by MAX_BACKFILL_WINDOWS) and log a WARNING
    containing 'backfill exhausted', rather than hang. Here EBI Search
    re-discovers the SAME already-failed accession on every window (no new
    accession ever appears), so no genuine progress is possible any more —
    that is functionally equivalent to exhaustion, and must NOT raise (only
    a cap-hit where EBI is still surfacing new, never-tried accessions
    raises — see test_cap_hit_with_entries_still_available_raises_runtime_error)."""

    def always_fails(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXDNEVERWORKS"}]}
        raise RuntimeError("permanent failure")

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(always_fails))
    with caplog.at_level(logging.WARNING):
        records = ingester._fetch_page('"chronic pain"', 0, 50)

    assert records == []
    assert any("backfill exhausted" in rec.message for rec in caplog.records)


def test_samples_count_zero_is_preserved_not_treated_as_missing(ingester, monkeypatch):
    """samplesCount: 0 is a genuine count and must not collapse to -1 the way
    samplesCount: None (unknown) does."""

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXDZERO"}, {"id": "PXDNONE"}]}
        acc = url.rstrip("/").split("/")[-1]
        samples = {"PXDZERO": 0, "PXDNONE": None}[acc]
        return {
            "accession": acc, "title": acc, "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2021-01-01",
            "samplesCount": samples, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))
    records = {r.accession: r for r in ingester._fetch_page('"q"', 0, 2)}
    assert records["PXDZERO"].n_samples == 0
    assert records["PXDNONE"].n_samples == -1


# ---------------------------------------------------------------------------
# Task 7 review follow-ups.
#
# Important 1: hitting MAX_BACKFILL_WINDOWS while EBI Search still has
#              entries must RAISE, not return a short page — a short page is
#              indistinguishable from genuine exhaustion to search_all().
# Important 2: search_all must be safely re-callable on the same instance.
# Minor 3:     a persistently-failing accession must cost at most one
#              metadata-fetch attempt (2 HTTP calls: try + retry) per
#              _fetch_page invocation, not one per backfill window.
# Minor 4:     that per-call failure tracking must NOT survive across calls,
#              so a transient failure can still recover on a later call.
# ---------------------------------------------------------------------------


def test_cap_hit_with_entries_still_available_raises_runtime_error(ingester, monkeypatch):
    """The backfill cap must not silently truncate a run: if EBI Search keeps
    discovering entries (never signals genuine exhaustion) but every one of
    them fails PRIDE metadata, hitting MAX_BACKFILL_WINDOWS must RAISE with
    'backfill exhausted' in the message, not return a short page."""
    calls = {"n": 0}

    def always_fails(url, params=None, timeout=30):
        if "ebisearch" in url:
            calls["n"] += 1
            # A fresh accession every window: EBI Search is never actually
            # exhausted, metadata just keeps failing.
            return {"entries": [{"id": f"PXDFAIL{calls['n']}"}]}
        raise RuntimeError("404")

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(always_fails))
    with pytest.raises(RuntimeError, match="backfill exhausted"):
        ingester._fetch_page('"chronic pain"', 40000, 50)


def test_genuine_exhaustion_with_entries_never_filling_page_does_not_raise(ingester, monkeypatch):
    """A source that is genuinely exhausted (EBI Search returns no more
    entries) must still return whatever records it found and must NOT raise
    — only a cap-hit-while-entries-remain condition raises."""
    corpus = [f"PXD4{i:06d}" for i in range(30)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, set()))
    )
    records = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(records) == 30


def test_search_all_is_safely_re_callable_on_the_same_instance(ingester, monkeypatch):
    """A second search_all() on the same instance must return the same
    record count as the first, not trip base.py's zero-recall guard because
    every accession was already in _seen_accessions from the first call."""
    corpus = [f"PXD5{i:06d}" for i in range(30)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, set()))
    )
    first = ingester.search_all(["chronic pain"], [], {}, page_size=100)
    second = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(first) == 30
    assert len(second) == 30


def test_failed_accession_is_not_retried_every_backfill_window(ingester, monkeypatch):
    """A permanently-failing accession must be attempted at most once (2 HTTP
    calls: the attempt and its retry) per _fetch_page invocation, not once
    per backfill window (5 windows x 2 retries = 10+ wasted calls)."""
    calls = {"PXDBAD": 0}

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            start = params["start"]
            # PXDBAD is rediscovered by EBI Search every window (it was never
            # successfully emitted, so it is never marked seen); exactly one
            # new good accession per window keeps the page filling slowly
            # enough to force several backfill windows.
            return {"entries": [{"id": "PXDBAD"}, {"id": f"PXDGOOD{start}"}]}
        acc = url.rstrip("/").split("/")[-1]
        if acc == "PXDBAD":
            calls["PXDBAD"] += 1
            raise RuntimeError("permanent 500")
        return {
            "accession": acc, "title": acc, "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2020-01-01",
            "samplesCount": 1, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))
    records = ingester._fetch_page('"q"', 0, 5)

    assert len(records) == 5
    assert "PXDBAD" not in {r.accession for r in records}
    assert calls["PXDBAD"] == 2, (
        "a permanently-failing accession must cost exactly 2 metadata calls "
        "(1 attempt + 1 retry) for the whole _fetch_page invocation, not 2 "
        "per backfill window"
    )


def test_cross_call_transient_recovery_after_failure_in_a_prior_fetch_page_call(
    ingester, monkeypatch
):
    """Minor 3's per-call failure tracking must be call-scoped: an accession
    that fails during one _fetch_page call must still be eligible for retry
    on a later _fetch_page call, in case the failure was transient."""
    state = {"attempts": 0}

    def flaky_across_calls(url, params=None, timeout=30):
        if "ebisearch" in url:
            return {"entries": [{"id": "PXDTRANSIENT"}]}
        state["attempts"] += 1
        # Fails both attempts (try + retry) within the FIRST _fetch_page
        # call, then succeeds once requested again on a later call.
        if state["attempts"] <= 2:
            raise RuntimeError("transient outage")
        return {
            "accession": "PXDTRANSIENT", "title": "Recovered on a later call",
            "projectDescription": "", "organisms": [], "instruments": [],
            "publicationDate": "2023-01-01", "samplesCount": 2, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(flaky_across_calls))

    first = ingester._fetch_page('"q"', 0, 1)
    assert first == [], "still down during the first _fetch_page call"

    second = ingester._fetch_page('"q"', 0, 1)
    assert [r.accession for r in second] == ["PXDTRANSIENT"], (
        "a failure local to one _fetch_page call must not block a later call "
        "from picking the accession back up once it recovers"
    )


# ---------------------------------------------------------------------------
# Task 7 SECOND review round.
#
# Important finding 1: an EBI Search DISCOVERY call failure (the `_get` to
#     EBI_SEARCH_URL itself throwing) was being treated as silent, genuine
#     exhaustion (`except Exception: logger.error(...); break`), returning a
#     short page with NO RuntimeError. That bypasses all the cap/raise
#     machinery added for metadata failures and contradicts the module
#     docstring's invariant that a short page always means "no more
#     studies". The fix retries the discovery call once (mirroring
#     `_fetch_project_metadata`) and raises RuntimeError (message contains
#     "discovery failed") if it still fails, rather than returning short.
#
# Important finding 2: `_fetch_page` advanced its internal cursor by the
#     full EBI Search window size even when the inner resolution loop
#     stopped early after filling the page, while search_all advances its
#     own offset by len(batch) (the emitted record count). The next
#     `_fetch_page` call therefore started behind where discovery had
#     actually gotten to and re-requested already-consumed EBI Search
#     windows — roughly doubling EBI traffic on any page that backfills.
#     Fixed by persisting the discovery cursor as instance state
#     (`self._cursor`) across `_fetch_page` calls, resolving every
#     discovered accession in a window (no early stop), and starting
#     discovery at `max(offset, self._cursor)`.
# ---------------------------------------------------------------------------


def test_discovery_failure_mid_backfill_raises_runtime_error(ingester, monkeypatch):
    """(a) A discovery-call failure with work still pending must raise
    RuntimeError containing 'discovery failed', not return a short page.
    The discovery call must actually be retried once before giving up."""
    calls = {"ebi": 0}

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            calls["ebi"] += 1
            start = params["start"]
            if start == 0:
                # First window: 3 accessions, all resolve fine, but that is
                # short of the requested page_size=10, so a backfill window
                # is required next.
                return {"entries": [{"id": "PXDA"}, {"id": "PXDB"}, {"id": "PXDC"}]}
            # Every backfill window's discovery call fails outright — this
            # is the reviewer's exact reproduction (EBI Search HTTP request
            # throwing on the backfill window).
            raise RuntimeError("EBI Search 500")
        acc = url.rstrip("/").split("/")[-1]
        return {
            "accession": acc, "title": acc, "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2020-01-01",
            "samplesCount": 1, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))

    with pytest.raises(RuntimeError, match="discovery failed"):
        ingester._fetch_page('"chronic pain"', 0, 10)

    # window1 (1 call) + window2 attempt + window2 retry (2 calls) = 3.
    # The old code broke after a single failed attempt (2 EBI calls total,
    # no retry) and returned a short page with no exception at all.
    assert calls["ebi"] == 3, "discovery must be retried once before raising"


def test_discovery_success_with_zero_entries_is_genuine_exhaustion_no_raise(
    ingester, monkeypatch
):
    """(b) A discovery call that SUCCEEDS but reports zero entries is
    genuine exhaustion: return whatever records were built, do not raise."""

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            start = params["start"]
            if start == 0:
                return {"entries": [{"id": "PXDONLY"}]}
            return {"entries": []}  # genuine exhaustion, not a failure
        return {
            "accession": "PXDONLY", "title": "Only one", "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2020-01-01",
            "samplesCount": 1, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))

    records = ingester._fetch_page('"chronic pain"', 0, 10)
    assert [r.accession for r in records] == ["PXDONLY"]


def test_discovery_fails_once_then_recovers_on_retry_no_raise(ingester, monkeypatch):
    """(c) A discovery call that fails on the first attempt but succeeds on
    the retry must not raise and must return the full page."""
    calls = {"n": 0}

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            calls["n"] += 1
            if calls["n"] == 1:
                raise RuntimeError("transient EBI 503")
            return {"entries": [{"id": "PXDOK1"}, {"id": "PXDOK2"}]}
        acc = url.rstrip("/").split("/")[-1]
        return {
            "accession": acc, "title": acc, "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2020-01-01",
            "samplesCount": 1, "keywords": [],
        }

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(fake))

    records = ingester._fetch_page('"chronic pain"', 0, 2)
    assert {r.accession for r in records} == {"PXDOK1", "PXDOK2"}
    assert calls["n"] == 2, "must retry exactly once before succeeding"


def test_redundant_discovery_calls_eliminated_across_a_250_study_corpus(
    ingester, monkeypatch
):
    """(d) Regression for the cursor-drift bug: over a full search_all() run
    on the 250-accession/1-permanently-failing corpus, the EBI Search
    `start=` values requested must be strictly non-decreasing with no
    duplicates, and the total call count must be <= 4 (previously 6, with
    duplicated start values [0, 100, 100, 200, 200, 250])."""
    corpus = [f"PXD1{i:06d}" for i in range(250)]
    broken = {corpus[41]}
    starts_seen: list[int] = []

    base_fake = _paginated_corpus_get(corpus, broken)

    def instrumented(url, params=None, timeout=30):
        if "ebisearch" in url:
            starts_seen.append(params["start"])
        return base_fake(url, params, timeout)

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(instrumented))

    records = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(records) == 249
    assert corpus[41] not in {r.accession for r in records}
    assert len({r.accession for r in records}) == 249, "no duplicate records"

    assert starts_seen == sorted(starts_seen), "start values must be non-decreasing"
    assert len(starts_seen) == len(set(starts_seen)), (
        f"no start value may be requested twice, got {starts_seen}"
    )
    assert len(starts_seen) <= 4, (
        f"expected <= 4 EBI Search calls (previously 6), got {len(starts_seen)}: "
        f"{starts_seen}"
    )


def test_search_all_twice_on_same_instance_after_cursor_fix_returns_same_count(
    ingester, monkeypatch
):
    """(e) Instance reuse must still work now that the discovery cursor is
    persisted instance state: two search_all() calls on the same instance
    must each return the full corpus, not be starved by a cursor left over
    from the previous call."""
    corpus = [f"PXD6{i:06d}" for i in range(75)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, set()))
    )

    first = ingester.search_all(["chronic pain"], [], {}, page_size=100)
    second = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(first) == 75
    assert len(second) == 75


# ---------------------------------------------------------------------------
# Task 7 THIRD review round.
#
# Important finding: `_fetch_page` can return MORE records than `page_size`.
# Because the per-accession early `break` was removed from
# `_resolve_accessions` (fixing the earlier skip bug), a window resolved in
# full during backfill can produce more records than the caller asked for.
# `search_all` never trims, so `hard_limit` is silently overshot. The fix
# buffers the overshoot in `self._pending` rather than trimming (trimming
# would strand accessions already added to `_seen_accessions`), and drains
# it at the top of the next `_fetch_page` call.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("hard_limit", "expected"),
    [(5, 5), (10, 10), (None, 49)],
)
def test_hard_limit_is_respected_exactly(ingester, monkeypatch, hard_limit, expected):
    """(a) A 50-accession corpus with one failing accession at index 2 must
    yield exactly `hard_limit` records when set, and the full 49 recoverable
    records when `hard_limit` is None. Reproduces the reviewer's overshoot:
    hard_limit=5 -> 9 records, hard_limit=10 -> 19 records."""
    corpus = [f"PXD7{i:06d}" for i in range(50)]
    broken = {corpus[2]}

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )
    records = ingester.search_all(
        ["chronic pain"], [], {}, hard_limit=hard_limit, page_size=100
    )

    assert len(records) == expected


def test_no_record_lost_to_buffering_across_full_search_all(ingester, monkeypatch):
    """(b) With page_size=10 and scattered metadata failures across a
    50-accession corpus, a full search_all() run must recover every
    successfully-fetchable accession exactly once — buffering overshoot
    into `_pending` must not drop or duplicate anything."""
    corpus = [f"PXDA{i:06d}" for i in range(50)]
    broken = {corpus[3], corpus[17], corpus[33], corpus[48]}

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )
    records = ingester.search_all(["chronic pain"], [], {}, page_size=10)

    accessions = [r.accession for r in records]
    assert set(accessions) == set(corpus) - broken
    assert len(accessions) == len(set(accessions)), "no duplicate records"
    assert len(records) == 46


def test_fetch_page_never_returns_more_than_page_size(ingester, monkeypatch):
    """(c) Direct check on the exact pattern that used to overshoot
    (hard_limit=5 -> 9 records): `_fetch_page` itself must never hand back
    more records than the `page_size` it was called with."""
    corpus = [f"PXDC{i:06d}" for i in range(50)]
    broken = {corpus[2]}

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )
    batch = ingester._fetch_page('"chronic pain"', 0, 5)

    assert len(batch) <= 5


def test_fetch_page_returns_full_page_while_pending_has_buffered_records(
    ingester, monkeypatch
):
    """(d) While `_pending` holds buffered overshoot, `_fetch_page` must
    return a FULL page — never a short one — so search_all keeps
    paginating instead of misreading buffered work as exhaustion. Only once
    pending is drained and discovery is genuinely exhausted may a page fall
    short of `page_size`."""
    corpus = [f"PXD9{i:06d}" for i in range(15)]
    broken = {corpus[2]}

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )

    all_records: list = []

    first = ingester._fetch_page('"chronic pain"', 0, 5)
    assert len(first) == 5
    assert ingester._pending, "the fix must buffer the overshoot, not discard it"
    all_records.extend(first)

    # Pending is non-empty: this call must return a FULL page, not a short
    # one, even though it is served partly from the buffer.
    second = ingester._fetch_page('"chronic pain"', len(all_records), 5)
    assert len(second) == 5
    all_records.extend(second)

    third = ingester._fetch_page('"chronic pain"', len(all_records), 5)
    all_records.extend(third)

    assert len(all_records) == 14, "genuine exhaustion must still yield every recoverable record"
    assert corpus[2] not in {r.accession for r in all_records}


def test_instance_reuse_resets_pending_buffer(ingester, monkeypatch):
    """(e) A prior run that stopped early via `hard_limit` (leaving
    `_pending` non-empty, since the overshoot beyond the cap is never
    drained) must not leak that buffered overshoot into a later
    `search_all()` call on the same instance."""
    corpus = [f"PXDD{i:06d}" for i in range(50)]
    broken = {corpus[2]}

    monkeypatch.setattr(
        ProteomicsIngester, "_get", staticmethod(_paginated_corpus_get(corpus, broken))
    )

    first = ingester.search_all(["chronic pain"], [], {}, hard_limit=5, page_size=100)
    assert len(first) == 5
    assert ingester._pending, "hard_limit must stop before the buffered overshoot is drained"

    second = ingester.search_all(["chronic pain"], [], {}, page_size=100)
    assert len(second) == 49
    assert ingester._pending == [], "reuse must reset the pending buffer"


# ---------------------------------------------------------------------------
# Task 7 FOURTH review round.
#
# Important finding: `_fetch_page` passed the FULL `page_size` (not the
# reduced remaining need) as `want` to `_resolve_new_page`, whose single
# `want` parameter drove both the stop condition/cap-raise gate AND the EBI
# Search discovery window size. When `_pending` already banked some records
# (normal after any overshoot), `_resolve_new_page` chased a full page of NEW
# records it did not actually need, burning through MAX_BACKFILL_WINDOWS
# under a degraded metadata-failure rate and raising "backfill exhausted"
# even though only a couple more records were required. The fix decouples
# the target (`want`, credited for pending work already banked) from the
# discovery window size (stays `page_size`-aligned for grid stability).
# ---------------------------------------------------------------------------


def _every_nth_succeeds_corpus_get(corpus: list[str], n: int):
    """Build a fake `_get` serving paginated EBI Search windows, where only
    every `n`-th accession's PRIDE metadata fetch succeeds (1-indexed
    position within `corpus`); all others permanently fail."""

    def fake(url, params=None, timeout=30):
        if "ebisearch" in url:
            start = params["start"]
            size = params["size"]
            chunk = corpus[start:start + size]
            return {"entries": [{"id": a} for a in chunk]}
        acc = url.rstrip("/").split("/")[-1]
        idx = corpus.index(acc)
        if (idx + 1) % n != 0:
            raise RuntimeError("404")
        return {
            "accession": acc, "title": f"Study {acc}", "projectDescription": "",
            "organisms": [], "instruments": [], "publicationDate": "2022-01-01",
            "samplesCount": 4, "keywords": [],
        }

    return fake


def test_pending_credit_prevents_spurious_backfill_exhausted_raise(ingester, monkeypatch):
    """(a) THE REGRESSION: with 8 records already banked in `_pending`,
    `_fetch_page("q", 8, 10)` only needs 2 NEW records. Against a corpus
    where just every 10th accession's metadata succeeds, 2 new records are
    reachable well within MAX_BACKFILL_WINDOWS extra windows — but the old
    code demanded a FULL 10 new records (ignoring the 8 already banked) and
    raised 'backfill exhausted'. Must return exactly 10 records, no raise."""
    corpus = [f"PXDN{i:06d}" for i in range(500)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get",
        staticmethod(_every_nth_succeeds_corpus_get(corpus, 10)),
    )

    ingester._pending = [
        ingester._to_record({
            "accession": f"PXDPEND{i}", "title": f"Pending {i}",
            "projectDescription": "", "organisms": [], "instruments": [],
            "publicationDate": "2020-01-01", "samplesCount": 1, "keywords": [],
        })
        for i in range(8)
    ]

    records = ingester._fetch_page('"q"', 8, 10)

    assert len(records) == 10


def test_empty_pending_genuine_cap_hit_still_raises(ingester, monkeypatch):
    """(b) Pending credit must not mask a genuine cap-hit: with `_pending`
    EMPTY and the same every-10th-succeeds corpus, a full page_size=10 worth
    of NEW records is NOT reachable within MAX_BACKFILL_WINDOWS extra
    windows, so `_fetch_page` must still raise 'backfill exhausted'."""
    corpus = [f"PXDM{i:06d}" for i in range(500)]

    monkeypatch.setattr(
        ProteomicsIngester, "_get",
        staticmethod(_every_nth_succeeds_corpus_get(corpus, 10)),
    )

    with pytest.raises(RuntimeError, match="backfill exhausted"):
        ingester._fetch_page('"q"', 0, 10)


def test_discovery_window_size_stays_page_size_not_reduced_need(ingester, monkeypatch):
    """(c) The EBI Search discovery window size must stay `page_size`
    (grid-stable), even though the TARGET (`want`) is reduced by pending
    credit. Instrument `_get` to record every `size=` passed to EBI Search
    during a `search_all()` run that triggers pending (backfill overshoot),
    and assert every recorded size equals page_size."""
    corpus = [f"PXDS{i:06d}" for i in range(60)]
    broken = {corpus[2]}
    sizes_seen: list[int] = []

    base_fake = _paginated_corpus_get(corpus, broken)

    def instrumented(url, params=None, timeout=30):
        if "ebisearch" in url:
            sizes_seen.append(params["size"])
        return base_fake(url, params, timeout)

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(instrumented))

    records = ingester.search_all(["chronic pain"], [], {}, page_size=10)

    assert len(records) == 59
    assert sizes_seen, "EBI Search must have been called at least once"
    assert all(size == 10 for size in sizes_seen), (
        f"discovery window size must stay page_size=10 for grid stability, got {sizes_seen}"
    )


def test_250_corpus_one_drop_guarantee_still_holds_after_decoupling(ingester, monkeypatch):
    """(d) The 250-accession/1-permanent-drop guarantee must still hold
    exactly after decoupling `want` from `window_size`: 249 records, 0
    duplicates, EBI `start` values [0, 100, 200, 250]."""
    corpus = [f"PXD1{i:06d}" for i in range(250)]
    broken = {corpus[41]}
    starts_seen: list[int] = []

    base_fake = _paginated_corpus_get(corpus, broken)

    def instrumented(url, params=None, timeout=30):
        if "ebisearch" in url:
            starts_seen.append(params["start"])
        return base_fake(url, params, timeout)

    monkeypatch.setattr(ProteomicsIngester, "_get", staticmethod(instrumented))

    records = ingester.search_all(["chronic pain"], [], {}, page_size=100)

    assert len(records) == 249
    assert corpus[41] not in {r.accession for r in records}
    assert len({r.accession for r in records}) == 249, "no duplicate records"
    assert starts_seen == [0, 100, 200, 250]
