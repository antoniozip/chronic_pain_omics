"""Proteomics ingester — EBI Search for discovery, PRIDE Archive for metadata.

The PRIDE v2 /projects endpoint accepts a `keyword` parameter and ignores it:
querying with and without it returns byte-identical results, so a keyword search
silently degrades to "the most recently published PRIDE datasets, any topic".
Discovery therefore goes through the EBI Search API, which filters correctly.

EBI Search returns empty `title` and `organism_name` fields for the pride domain,
so metadata is fetched per accession from PRIDE.

A single transient PRIDE metadata failure must never silently truncate a run:
base.OmicsIngester.search_all() treats a page shorter than requested as "source
exhausted" and stops. Each failed metadata fetch is retried once; a persistent
failure drops that accession only from THIS window and is logged, but
`_fetch_page` then backfills with additional, real EBI Search windows until it
either fills `page_size` or the source is genuinely exhausted — so a page
short of `page_size` once again means "no more studies", not "one dropped
accession", which is the invariant search_all's exhaustion check relies on.

The discovery call itself (EBI Search) is just as capable of failing as a
PRIDE metadata fetch. If it throws, `_fetch_discovery_window` retries it once;
a failure that survives the retry is escalated to a RuntimeError (message
contains "discovery failed") rather than returned as a short page — returning
a short page here would be indistinguishable, to search_all's exhaustion
check, from genuine exhaustion, and would silently truncate the run exactly
like an unhandled metadata failure would. Only a discovery call that
*succeeds* and reports zero entries is genuine exhaustion.

`_fetch_page` bounds the backfill loop at `MAX_BACKFILL_WINDOWS` extra windows
so a pathologically broken source cannot hang it forever. But hitting that cap
is NOT always the same thing as genuine exhaustion: if EBI Search keeps
discovering accessions PRIDE metadata has never been tried against, and every
one of them fails, returning a short page would let search_all silently read
it as "no more studies" — exactly the truncation-by-log-line failure class
this module exists to eliminate. `_fetch_page` therefore raises RuntimeError
when the cap is hit and the most recent window still turned up at least one
accession that had not already been tried (successfully or not) earlier in
this same call. It returns a short (or empty) page without raising when EBI
Search itself ran out of entries, or when every accession the cap-hitting
windows kept re-surfacing had already been tried this call — both mean no
further progress is possible regardless of how many more windows are spent.

Because backfill can pull in accessions beyond the caller's requested window,
and search_all advances its own offset by len(batch), every accession this
ingester ever emits is tracked in `_seen_accessions` so it is never emitted
twice across successive `_fetch_page` calls. That tracker is instance state,
so `search_all` is overridden here to reset it on every call — otherwise a
second `search_all()` on the same instance would treat every accession as
already emitted and trip base.py's zero-recall guard.

The EBI Search discovery cursor is likewise instance state (`_cursor`),
persisted across `_fetch_page` calls and reset alongside `_seen_accessions`
in the `search_all` override. Every discovery window fetched by `_fetch_page`
is resolved in full — no accession is left un-examined mid-window — so the
cursor can always be advanced by exactly the window size with nothing lost.
Without persistence, each new `_fetch_page` call (driven by search_all's own
`offset`, which only counts *emitted* records) would re-request EBI Search
windows this ingester had already consumed while backfilling past a dropped
accession, roughly doubling EBI Search traffic on any page that backfills.

Resolving a window in full (rather than stopping once `records` reaches
`page_size`) means a window can produce MORE records than the caller asked
for. Those records must not be discarded — every one of them is already in
`_seen_accessions` and would be permanently lost if trimmed, silently
undershooting the caller's `hard_limit` on a LATER page instead of exactly
honoring it now. So the overshoot is buffered in `self._pending` (instance
state, reset alongside `_seen_accessions`/`_cursor`) and drained into the
front of the next `_fetch_page` call's result, in order, before any new
discovery happens. `_fetch_page` therefore never returns more than
`page_size` records, and base.py's `hard_limit` check in `search_all` (which
never trims) is never overshot. Because a non-empty `_pending` means real,
already-fetched records are waiting, `_fetch_page` always returns a FULL
page while `_pending` is non-empty — only once it is drained AND discovery
genuinely reports no more entries does a page fall short of `page_size`,
preserving search_all's exhaustion contract (`len(batch) < page_size` means
"no more studies").

Crediting `_pending` against the target must not also shrink the EBI Search
discovery window: `_fetch_page` calls `_resolve_new_page` with the remaining
need (`page_size - len(records)` after draining `_pending`) as the TARGET,
but keeps the discovery WINDOW SIZE pinned to the full `page_size`. Without
this split, a non-empty `_pending` (normal after any overshoot) made
`_resolve_new_page` demand a full page of NEW records it did not actually
need, burning through `MAX_BACKFILL_WINDOWS` under a degraded metadata-
failure rate and raising "backfill exhausted" even when only a couple more
records were required to complete the page. The window size stays
`page_size`-aligned so the `start=` grid `_resolve_new_page` walks never
drifts across calls, while the target correctly reflects only the work still
outstanding.
"""

from __future__ import annotations

import logging
from typing import Any

from . import register_ingester
from .base import DatasetRecord, OmicsIngester

logger = logging.getLogger(__name__)

EBI_SEARCH_URL = "https://www.ebi.ac.uk/ebisearch/ws/rest/pride"
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"
PRIDE_PROJECT_URL = f"{PRIDE_API_BASE}/projects"
PRIDE_STUDY_UI = "https://www.ebi.ac.uk/pride/archive/projects/"

# Safety cap on extra EBI Search windows _fetch_page will consume trying to
# backfill a page after drops. Bounds the loop even if every accession in the
# source keeps failing metadata (or keeps coming back already-seen).
MAX_BACKFILL_WINDOWS = 5


@register_ingester("proteomics")
class ProteomicsIngester(OmicsIngester):
    """Queries PRIDE Archive (via EBI Search) for chronic-pain proteomics datasets."""

    modality = "proteomics"

    def __init__(self) -> None:
        super().__init__()
        # Accessions already emitted by this ingester instance, across all
        # _fetch_page calls. Required so backfill windows (which read beyond
        # the caller's [offset, offset+page_size) window) never re-emit an
        # accession search_all already received.
        self._seen_accessions: set[str] = set()
        # Position in the EBI Search result stream this ingester has already
        # consumed, across all _fetch_page calls. Persisted (rather than
        # reset to the caller's `offset` on every call) so a page that
        # backfills past `offset + page_size` never causes the NEXT
        # `_fetch_page` call to re-request windows already fetched (see
        # module docstring).
        self._cursor: int = 0
        # Records resolved during a prior _fetch_page call that overshot the
        # page_size requested at the time (because window resolution never
        # stops early — see module docstring). Drained, in order, at the top
        # of the next _fetch_page call before any new discovery happens, so
        # no record already added to _seen_accessions is ever lost.
        self._pending: list[DatasetRecord] = []
        # Set once a discovery call has SUCCEEDED and reported zero entries
        # (genuine exhaustion, see `_fetch_discovery_window`). Buffering
        # overshoot into `_pending` can split what used to be one page's
        # worth of resolution across two `_fetch_page` calls; without this
        # flag, the second call would re-probe the same already-exhausted
        # EBI Search position, doubling discovery traffic right at the tail
        # of a run. Once true, no further discovery calls are attempted this
        # instance's lifetime (until `search_all` resets it).
        self._exhausted: bool = False
        # Total hits reported by EBI Search for this query, learned from the
        # first successful discovery response and used to stop paging before
        # the endpoint 400s past the end.
        self._hit_count: int | None = None

    def _build_query(self, terms: list[str], species_terms: list[str]) -> str:
        return " OR ".join(f'"{t}"' for t in terms)

    def search_all(self, *args: Any, **kwargs: Any) -> list[DatasetRecord]:
        """Reset the seen-accession tracker and discovery cursor, then delegate.

        Both are instance state (see module docstring) so that backfill
        windows within a single `search_all()` run never re-emit an
        accession, or re-request an EBI Search window, `_fetch_page` already
        consumed. Left un-reset, a second `search_all()` call on the same
        ingester instance would find every accession already "seen", emit
        nothing, and trip base.py's zero-recall guard with a misleading
        "returned 0 records" error — exactly the silent/false zero this
        module's docstring warns about. Resetting here makes the ingester
        safely re-callable without changing any other base-class behaviour.
        """
        self._seen_accessions = set()
        self._cursor = 0
        self._pending = []
        self._exhausted = False
        self._hit_count = None
        return super().search_all(*args, **kwargs)

    def _fetch_discovery_window(
        self, query: str, start: int, page_size: int
    ) -> list[str] | None:
        """Fetch one EBI Search discovery window, retrying once on failure.

        Mirrors `_fetch_project_metadata`'s retry-once policy so a transient
        EBI Search outage cannot masquerade as genuine exhaustion any more
        than a transient PRIDE metadata outage can.

        Returns:
            The window's accession ids (an empty list means EBI Search
            itself reports no more entries — genuine exhaustion). ``None``
            means the discovery call failed on both the initial attempt and
            the retry; the caller must not treat that the same as an empty
            list.
        """
        # EBI Search answers `start` beyond the last hit with HTTP 400, not with
        # an empty page, so the "a short page means exhaustion" rule below can
        # never fire on an exact multiple and the run dies one page past the
        # end. hitCount, which every response carries, is the honest bound: a
        # pain query returning exactly 43 hits must stop at 43 rather than ask
        # for start=43 and read the 400 as source failure.
        if self._hit_count is not None and start >= self._hit_count:
            return []

        for attempt in range(2):
            try:
                data = self._get(
                    EBI_SEARCH_URL,
                    params={
                        "query": query, "format": "json",
                        "size": page_size, "start": start,
                    },
                )
            except Exception as exc:
                if attempt == 0:
                    continue
                logger.error(
                    "[%s] EBI Search discovery failed at start=%d after retry: %s",
                    self.modality, start, exc,
                )
                return None
            hits = data.get("hitCount")
            if isinstance(hits, int) and self._hit_count is None:
                self._hit_count = hits
                logger.info("[%s] EBI Search reports %d hits", self.modality, hits)
            return [
                e["id"] for e in data.get("entries", []) or []
                if isinstance(e, dict) and e.get("id")
            ]
        return None

    def _resolve_accessions(
        self,
        accessions: list[str],
        records: list[DatasetRecord],
        attempted: set[str],
        failed: set[str],
    ) -> None:
        """Fetch PRIDE metadata for every accession in one discovery window.

        Every accession is examined — this never stops early once `records`
        reaches the caller's page size. Consuming the whole window on every
        fetch is what lets `_fetch_page` advance the persisted discovery
        cursor by exactly the window size with nothing left un-examined
        (see module docstring on `_cursor`); stopping early would either
        lose track of the un-examined tail or force it to be re-fetched.
        """
        for accession in accessions:
            if accession in self._seen_accessions or accession in failed:
                continue
            attempted.add(accession)
            meta = self._fetch_project_metadata(accession)
            if meta is None:
                failed.add(accession)
                continue
            records.append(self._to_record(meta))
            self._seen_accessions.add(accession)

    def _fetch_page(self, query: str, offset: int, page_size: int) -> list[DatasetRecord]:
        """Return at most `page_size` records, buffering any overshoot.

        Window resolution (`_resolve_accessions`) never stops early, so a
        single window can yield more records than the caller asked for.
        Trimming that overshoot would strand accessions already committed
        to `_seen_accessions`, so instead it is drained from — and buffered
        into — `self._pending` (see module docstring). Any records already
        pending from a previous call are emitted first, before any new
        discovery happens.

        A non-empty `self._pending` after this call means real, already-
        fetched work is waiting, so this always returns a FULL page in that
        case — never a short one — preserving search_all's exhaustion
        contract (`len(batch) < page_size` means genuine exhaustion).

        `_resolve_new_page` is called with the REMAINING need
        (`page_size - len(records)`, i.e. `page_size` minus whatever
        `_pending` already contributed) as its target, but with the discovery
        WINDOW SIZE held at the full `page_size` regardless. These two used
        to be the same parameter, which meant a non-empty `_pending` made
        `_resolve_new_page` chase a full page of NEW records it did not
        actually need, burning through `MAX_BACKFILL_WINDOWS` under a
        degraded metadata-failure rate and raising "backfill exhausted" even
        when only a couple more records were required. Decoupling them fixes
        that: the target now credits work already banked in `_pending`, while
        the window size stays `page_size`-aligned so the EBI Search `start=`
        grid never drifts (see `_resolve_new_page`'s docstring and the test
        asserting the start sequence `[0,100,200,250]` with no duplicates).
        """
        records: list[DatasetRecord] = []
        if self._pending:
            take = min(len(self._pending), page_size)
            records.extend(self._pending[:take])
            self._pending = self._pending[take:]
            if len(records) >= page_size:
                return records[:page_size]

        need = page_size - len(records)
        new_records = self._resolve_new_page(query, offset, need, page_size)
        records.extend(new_records)

        if len(records) > page_size:
            self._pending.extend(records[page_size:])
            records = records[:page_size]

        return records

    def _resolve_new_page(
        self, query: str, offset: int, want: int, window_size: int
    ) -> list[DatasetRecord]:
        """Discover and resolve fresh accessions until at least `want` records
        are built (or the source is exhausted/capped).

        `want` and `window_size` are deliberately separate parameters. `want`
        drives the `while len(records) < want` stop condition and the
        cap-hit raise gate (including its message) below — it is the caller's
        actual remaining need, already credited for any work `_fetch_page`
        banked in `_pending`. `window_size` drives only the EBI Search
        discovery call (`_fetch_discovery_window(query, window_start,
        window_size)`) and stays pinned to the caller's `page_size` so the
        `start=` grid never drifts across calls, regardless of how small
        `want` is. The cursor advance (`self._cursor = window_start +
        len(accessions)`) is unaffected by either: it always advances by the
        actual accession count a window returned.

        May return MORE than `want` records: a discovery window is always
        resolved in full (see `_resolve_accessions`), so the window that
        finally reaches `want` can overshoot it. The caller (`_fetch_page`)
        is responsible for buffering any overshoot; this method's contract
        is only "at least `want`, unless genuinely exhausted or capped".

        If `self._exhausted` is already set (a prior call's discovery
        succeeded with zero entries — genuine exhaustion), returns
        immediately with no discovery call at all: buffering overshoot into
        `_pending` can split what used to be one page's worth of resolution
        across two `_fetch_page` calls, and without this short-circuit the
        second call would re-probe the exact EBI Search position already
        known empty, doubling discovery traffic at the tail of a run.
        """
        if self._exhausted:
            return []

        records: list[DatasetRecord] = []
        attempted: set[str] = set()
        # Call-scoped only: an accession that fails metadata fetch within
        # THIS call is skipped for the rest of this call (so a persistent
        # failure isn't retried on every backfill window), but a later
        # _fetch_page/search_all call gets a fresh `failed` set and will try
        # it again in case the failure was transient.
        failed: set[str] = set()
        extra_windows = 0
        first_window = True
        # Whether the most recently fetched window turned up at least one
        # accession not already tried (successfully or not) earlier in this
        # call. Only that signals genuine, ongoing work remains — a window
        # that only re-surfaces already-seen/already-failed accessions means
        # no further progress is possible no matter how many more windows
        # are spent, which is functionally exhaustion even if EBI Search
        # technically returned a non-empty entries list.
        last_window_made_progress = False

        # max() guards a caller-supplied offset ahead of the persisted
        # cursor (never happens via search_all, but keeps this correct for
        # any other caller); in the normal case self._cursor already leads
        # offset because dropped accessions make the cursor advance faster
        # than emitted records.
        window_start = max(offset, self._cursor)

        while len(records) < want:
            if not first_window:
                extra_windows += 1
                if extra_windows > MAX_BACKFILL_WINDOWS:
                    logger.warning(
                        "[%s] backfill exhausted after %d extra EBI Search windows; "
                        "returning %d of %d requested records",
                        self.modality, MAX_BACKFILL_WINDOWS, len(records), want,
                    )
                    if last_window_made_progress:
                        # EBI Search is still discovering accessions that have
                        # never been tried this call. A short page here is
                        # indistinguishable to search_all() from "no more
                        # studies", so it would truncate the run silently.
                        # Escalate instead of returning a short page.
                        raise RuntimeError(
                            f"[{self.modality}] backfill exhausted after "
                            f"{MAX_BACKFILL_WINDOWS} extra EBI Search windows at "
                            f"offset {offset}: only {len(records)} of {want} "
                            f"requested records were built from {len(attempted)} "
                            "attempted accessions, and EBI Search is still "
                            "discovering new accessions. PRIDE metadata fetches "
                            "are persistently failing; this is not genuine "
                            "source exhaustion."
                        )
                    break
            first_window = False

            accessions = self._fetch_discovery_window(query, window_start, window_size)
            if accessions is None:
                # The discovery call itself failed twice (see
                # _fetch_discovery_window). We only ever get here while
                # len(records) < want, i.e. more work is always pending, so
                # a short page would be indistinguishable from genuine
                # exhaustion to search_all(). Escalate instead.
                raise RuntimeError(
                    f"[{self.modality}] EBI Search discovery failed at "
                    f"offset={offset}, cursor={window_start} after retry: "
                    f"{len(records)} of {want} requested records were "
                    "built before the failure. This is not genuine source "
                    "exhaustion."
                )

            if not accessions:
                # Genuine exhaustion: no more accessions at all. Remember
                # this so a later _fetch_page call (e.g. one drained mostly
                # from _pending) never re-probes this same EBI Search
                # position, which it would otherwise have no way to know is
                # already known to be empty.
                self._exhausted = True
                break

            attempted_before = len(attempted)
            self._resolve_accessions(accessions, records, attempted, failed)

            # Persist the cursor past this whole window: _resolve_accessions
            # never stops early, so nothing from `accessions` is left
            # un-examined for a later call to lose track of.
            self._cursor = window_start + len(accessions)
            window_start = self._cursor

            last_window_made_progress = len(attempted) > attempted_before

        if len(records) < want and len(attempted) > len(records):
            logger.warning(
                "[%s] page shrank: EBI Search discovered %d accessions but only %d "
                "PRIDE metadata records were built; some accessions were dropped",
                self.modality, len(attempted), len(records),
            )

        return records

    def _fetch_project_metadata(self, accession: str) -> dict[str, Any] | None:
        """Fetch full project metadata; EBI Search does not populate title/organism.

        Retries once on failure before giving up, so a single transient error
        cannot masquerade as a genuinely short page (see module docstring).
        """
        for attempt in range(2):
            try:
                return self._get(f"{PRIDE_PROJECT_URL}/{accession}")
            except Exception as exc:
                if attempt == 0:
                    continue
                logger.warning(
                    "metadata dropped for accession %s after retry: %s", accession, exc
                )
                return None
        return None

    def _to_record(self, item: dict[str, Any]) -> DatasetRecord:
        accession = item.get("accession", "")

        species = [
            o.get("name", "") for o in item.get("organisms", []) or [] if isinstance(o, dict)
        ]
        platform = "; ".join(
            i.get("name", "") for i in item.get("instruments", []) or [] if isinstance(i, dict)
        )

        pub_date = item.get("publicationDate", "") or item.get("submissionDate", "")

        # samplesCount is frequently None in PRIDE v2; -1 means "unknown", not
        # "zero" — a genuine samplesCount: 0 must be preserved as 0.
        sc = item.get("samplesCount")
        sample_count = sc if sc is not None else -1

        keywords = item.get("keywords", []) or []
        if isinstance(keywords, str):
            keywords = [keywords]

        return DatasetRecord(
            modality=self.modality,
            source="PRIDE",
            accession=accession,
            title=item.get("title", ""),
            description=item.get("projectDescription", ""),
            species=species,
            n_samples=int(sample_count),
            platform=platform or "Mass spectrometry",
            year=str(pub_date)[:4],
            url=f"{PRIDE_STUDY_UI}{accession}",
            keywords=keywords,
        )
