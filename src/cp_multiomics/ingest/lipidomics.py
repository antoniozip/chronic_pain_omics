"""Lipidomics ingester — the lipidomics subset of MetaboLights.

MetaboLights hosts metabolomics and lipidomics together, so this is
`MetabolomicsIngester` with its modality filter inverted: same EBI Search
discovery, same parsing, same `MTBLS`-only rule, opposite verdict on the
lipidomics keywords. The class docstring always claimed to "reuse
MetabolomicsIngester's raw fetch logic"; before 2026-08-15 it duplicated that
logic instead, and duplicated its two defects with it — a query parameter
MetaboLights ignores, and a response shape the endpoint does not return. See
`ingest/metabolomics.py` for both.

Inverting one predicate, rather than keeping two keyword lists in step by hand,
is also what guarantees every study is claimed by exactly one of the two
ingesters: no study can be dropped by both or counted by both.
"""

from __future__ import annotations

import logging

from . import register_ingester
from .metabolomics import MetabolomicsIngester

logger = logging.getLogger(__name__)

@register_ingester("lipidomics")
class LipidomicsIngester(MetabolomicsIngester):
    """Queries MetaboLights for lipidomics-flagged datasets on pain.

    `_exclude_lipid_keywords` is deliberately *not* overridden: inheriting the
    parent's set is what makes the two predicates exact complements. A private
    keyword list here could drift from the parent's, and a study matching one
    list but not the other would be claimed by both ingesters or by neither.
    """

    modality = "lipidomics"

    def _keep(self, text: str) -> bool:
        """Keep exactly what MetabolomicsIngester discards."""
        return self._is_lipidomics(text)
