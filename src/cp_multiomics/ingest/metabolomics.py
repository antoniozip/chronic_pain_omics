"""Metabolomics ingester — MetaboLights via EBI Search, plus Metabolomics Workbench.

**Two repositories, one manifest.** Until 2026-08-27 this ingester searched
MetaboLights alone, and the metabolomic arm reached three analysable studies
with one of them carrying 90% of the pooled weight. Metabolomics Workbench is
the field's other major deposition site and was never queried. Both now feed
`data/raw/metabolomics/datasets.jsonl`, distinguished by `source`, because
`04_harmonize.py` reads exactly that path per modality and a systematic
review's corpus is the union of its sources. See `metabolomics_workbench.py`
for why discovery there enumerates rather than queries.

## MetaboLights

The MetaboLights `/ws/studies` endpoint accepts a `query` parameter and ignores
it. Verified on 2026-08-15: the responses to no query, `query=pain` and
`query=zzzznonsense` are byte-identical (md5 e961b3d6...), and all three are the
complete list of 3,281 accessions. It also returns that list as bare strings,
not study objects, so the previous implementation's `isinstance(s, dict)` guard
skipped every element and the ingester emitted nothing — which is why
`data/raw/metabolomics/datasets.jsonl` sat at 0 bytes and the arm was recorded
as blocked on a "broken keyword API" needing manual curation.

Both defects are the PRIDE `/projects` story again (see `ingest/proteomics.py`),
and the fix is the same: discovery goes through EBI Search, which filters
correctly. Unlike the pride domain, the metabolights domain returns full
metadata inline, so there is no second per-accession metadata call and none of
proteomics' backfill machinery is needed here.

**The domain mixes studies and compounds.** It indexes 36,416 entries against
3,281 studies: `MTBLS*` are studies, `MTBLC*` are reference compounds. A pain
query returns both — `MTBLC27808` is heroin, `MTBLC3216` is buprenorphine —
and emitting a compound as a dataset would put a row through the pipeline that
has no samples and no design. Only `MTBLS` accessions are kept.

**Retrieval is deliberately broader than the configured term list.** The four
configured phrases ("chronic pain", "neuropathic pain", ...) return 5 entries,
3 of them studies, because MetaboLights titles say "Recurrent Pelvic Pain" or
"Chronic Craniofacial Pain" rather than the exact bigrams. Requiring those
bigrams is an artefact of query syntax, not an inclusion criterion. A bare
`pain` anchor is therefore ORed in, giving 41 studies to screen. That is the
PRISMA order — retrieve with high sensitivity, exclude at triage with a reason
code — and it is affordable here in a way it would not be for PRIDE, because
the whole repository is 3,281 studies.
"""

from __future__ import annotations

import logging
from typing import Any

from . import metabolomics_workbench as workbench
from . import register_ingester
from .base import DatasetRecord, OmicsIngester

logger = logging.getLogger(__name__)

EBI_SEARCH_URL = "https://www.ebi.ac.uk/ebisearch/ws/rest/metabolights"
METABOLIGHTS_STUDY_UI = "https://www.ebi.ac.uk/metabolights/"

# Fields the metabolights domain returns inline. Requesting a field the domain
# does not define yields an empty list rather than an error, so an entry with
# no organism is normal and not a failure.
EBI_FIELDS = (
    "name,description,organism,organism_part,technology_type,"
    "instrument_platform,study_design,study_factor,submission_date,publication_date"
)

# Study accessions, as opposed to MTBLC reference compounds. See module docstring.
STUDY_PREFIX = "MTBLS"


@register_ingester("metabolomics")
class MetabolomicsIngester(OmicsIngester):
    """Queries MetaboLights, via EBI Search, for metabolomics datasets on pain.

    Filters out lipidomics-specific studies; those are the LipidomicsIngester's.
    """

    modality = "metabolomics"
    # "lipidomic" (no s) is a prefix of "lipidomics", so this one entry covers
    # both spellings. LipidomicsIngester inherits this exact set and inverts
    # the verdict; if it kept its own list the two could disagree and a study
    # would be claimed by both ingesters or by neither.
    _exclude_lipid_keywords = {"lipidomic", "lipidome", "lipid profiling"}

    # Broad anchor ORed into the configured terms; see module docstring for why
    # the configured phrases alone under-retrieve on this repository.
    _broad_terms = ("pain",)

    def __init__(self) -> None:
        # hitCount from the first response. EBI Search answers a `start` beyond
        # the last hit with HTTP 400 rather than an empty page, so "a short page
        # means exhaustion" can never fire on an exact multiple and the run
        # would die one page past the end. This is the same defect fixed for
        # proteomics in d0ebf5e.
        self._hit_count: int | None = None
        # Our own position in the EBI result set. It cannot be search_all's
        # `offset`, which counts records *emitted*: this ingester drops
        # compounds and lipidomics studies, so a window of 100 entries may
        # yield 70 records, and reusing the emitted count as the next `start`
        # would re-read 30 entries and never reach the tail. Harmless while the
        # whole result set fits in one page, wrong the moment it does not.
        self._cursor = 0
        # Records parsed but not yet emitted. A window can yield more records
        # than the caller asked for, and the cursor has already moved past
        # them, so dropping the surplus would lose those studies outright.
        self._pending: list[DatasetRecord] = []
        # Metabolomics Workbench is drained once, after EBI Search exhausts.
        # Without the flag it would be re-enumerated on every subsequent page
        # — 4,507 studies and ~85 s per call — and the same studies emitted
        # repeatedly.
        self._workbench_done = False
        # Workbench exclusions, held for the PRISMA record.
        self._workbench_excluded: list[dict] = []

    def _build_query(self, terms: list[str], species_terms: list[str]) -> str:
        quoted = [f'"{t}"' for t in terms]
        return " OR ".join([*quoted, *self._broad_terms])

    def search_all(self, *args: Any, **kwargs: Any) -> list[DatasetRecord]:
        """Reset the hit-count bound, then delegate.

        Both are instance state, so a second call would otherwise inherit the
        first run's bound and cursor and stop early.
        """
        self._hit_count = None
        self._cursor = 0
        self._pending = []
        self._workbench_done = False
        self._workbench_excluded = []
        return super().search_all(*args, **kwargs)

    def _fetch_page(self, query: str, offset: int, page_size: int) -> list[DatasetRecord]:
        """Return up to `page_size` records, reading as many EBI windows as needed.

        `offset` is ignored in favour of `self._cursor`; see __init__. Windows
        are read until the page is full or EBI Search is exhausted, so a page
        falling short of `page_size` means "no more studies" — the invariant
        search_all's exhaustion check relies on — rather than "this window
        happened to be mostly compounds".
        """
        while len(self._pending) < page_size:
            if self._hit_count is not None and self._cursor >= self._hit_count:
                self._drain_workbench()
                break
            window = self._fetch_window(query, self._cursor, page_size)
            self._cursor += page_size
            if not window:
                self._drain_workbench()
                break
            self._pending.extend(self._parse_entries(window))

        page, self._pending = self._pending[:page_size], self._pending[page_size:]
        return page

    def _drain_workbench(self) -> None:
        """Append the Workbench records once EBI Search is exhausted.

        Both repositories feed one manifest, `data/raw/metabolomics/datasets.jsonl`,
        because `04_harmonize.py` reads exactly that path per modality and a
        systematic review's corpus is the union of its sources, distinguished
        by `source`. Draining here rather than in `search_all` keeps the
        "a short page means exhausted" invariant intact: the Workbench records
        enter `_pending` before that page is cut.

        Accession spaces do not collide (MTBLS* against ST*), so a study
        deposited in both repositories is caught at triage and recorded in
        `conf/analysis/superseded_studies.csv`, never dropped silently here.
        """
        if self._workbench_done:
            return
        self._workbench_done = True
        try:
            studies = workbench.fetch_all(self._get)
        except Exception as exc:
            # A failure here must not look like exhaustion of a repository
            # that was never reached: half a corpus silently becomes the
            # whole corpus, and the PRISMA count is wrong with no error.
            raise RuntimeError(
                f"[{self.modality}] Metabolomics Workbench discovery failed: {exc}"
            ) from exc
        kept, excluded = workbench.screen(studies, self.modality)
        self._workbench_excluded = excluded
        logger.info("[%s] Workbench: %d enumerated, %d kept, %d excluded",
                    self.modality, len(studies), len(kept), len(excluded))
        self._pending.extend(workbench.to_record(s, self.modality) for s in kept)

    def _fetch_window(self, query: str, start: int, page_size: int) -> list[dict]:
        """One EBI Search window, retried once. Empty means genuine exhaustion."""
        for attempt in range(2):
            try:
                data = self._get(
                    EBI_SEARCH_URL,
                    params={
                        "query": query, "format": "json",
                        "size": page_size, "start": start,
                        "fields": EBI_FIELDS,
                    },
                )
            except Exception as exc:
                if attempt == 0:
                    continue
                # Returning [] here would be indistinguishable, to
                # search_all's exhaustion check, from a genuinely empty page,
                # and would silently truncate the run.
                raise RuntimeError(
                    f"[{self.modality}] EBI Search discovery failed at "
                    f"start={start} after retry: {exc}"
                ) from exc

            hits = data.get("hitCount")
            if isinstance(hits, int) and self._hit_count is None:
                self._hit_count = hits
                logger.info("[%s] EBI Search reports %d hits", self.modality, hits)
            return data.get("entries", []) or []
        return []

    def _parse_entries(self, entries: list[dict]) -> list[DatasetRecord]:
        records: list[DatasetRecord] = []
        for e in entries:
            if not isinstance(e, dict):
                continue
            accession = e.get("id", "")
            if not accession.startswith(STUDY_PREFIX):
                continue        # reference compound, not a study

            f = e.get("fields", {}) or {}
            title = _first(f, "name")
            description = _first(f, "description")

            design = (f.get("study_design", []) or []) + (f.get("study_factor", []) or [])
            combined = " ".join([title, description, *design]).lower()
            if not self._keep(combined):
                continue

            # "reference compound" and "blank" appear as organism values on QC
            # entries; they are not species and must not reach harmonization.
            species = [
                o for o in dict.fromkeys(f.get("organism", []) or [])
                if o and o.lower() not in {"reference compound", "blank", "mixed sample"}
            ]

            platform = "; ".join(
                dict.fromkeys(
                    (f.get("technology_type", []) or []) + (f.get("instrument_platform", []) or [])
                )
            )
            year = (_first(f, "publication_date") or _first(f, "submission_date"))[:4]

            records.append(
                DatasetRecord(
                    modality=self.modality,
                    source="METABOLIGHTS",
                    accession=accession,
                    title=title,
                    description=description,
                    species=species,
                    # EBI Search carries no sample count. -1 means "unknown", which
                    # apply_filters passes through; 0 would be read as "fewer than
                    # min_samples" and silently drop every study.
                    n_samples=-1,
                    platform=platform or "Metabolomics",
                    year=year,
                    url=f"{METABOLIGHTS_STUDY_UI}{accession}",
                    keywords=list(dict.fromkeys(
                        (f.get("study_design", []) or []) + (f.get("study_factor", []) or [])
                    )),
                )
            )
        return records

    def _keep(self, text: str) -> bool:
        """Whether an entry belongs to this modality.

        MetaboLights hosts both; LipidomicsIngester inverts this rather than
        re-implementing the fetch, so a study is claimed by exactly one of the
        two ingesters and neither can silently swallow the other's studies.
        """
        return not self._is_lipidomics(text)

    def _is_lipidomics(self, text: str) -> bool:
        return any(kw in text for kw in self._exclude_lipid_keywords)


def _first(fields: dict, key: str) -> str:
    values = fields.get(key, []) or []
    return values[0] if values else ""
