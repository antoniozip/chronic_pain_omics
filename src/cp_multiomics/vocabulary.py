"""Single source of truth for the chronic-pain search vocabulary.

As of Stage 0 this vocabulary drives SCREENING (pipeline/02_screen.py auto-exclusion)
and the proteomics relevance scorer (scripts/analyze_proteomics.py). It is NOT yet wired
into retrieval: ingestion (pipeline/03_ingest_omics.py) and literature search
(pipeline/01_search_literature.py) still issue the original 4-term query. Widening
retrieval to the full vocabulary is deferred to Stage 1, where the resulting ~9x jump in
retrieved records is absorbed together with per-modality cohort expansion
(see plan/cohort_expansion_plan.md). Keeping the vocabulary in one place now means that
Stage-1 wiring is a config change, not a re-extraction.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

DEFAULT_VOCAB_PATH = (
    Path(__file__).resolve().parents[2] / "conf" / "search" / "pain_vocabulary.yaml"
)

# Short acronyms (<=4 chars) are prone to matching inside unrelated words
# (e.g. "CCI" inside "occipital"). For those we require a word boundary and
# allow an optional trailing "s" for the plural (e.g. "DRGs"), so we never
# silently drop a real study over a missing plural. Longer terms keep the
# original, more permissive substring behaviour. Patterns are compiled once
# per term and cached at module scope so the frozen dataclass never needs a
# mutable cache attribute.
SHORT_TERM_MAX_LEN = 4


@cache
def _short_term_pattern(term: str) -> re.Pattern[str]:
    """Compile (and cache) the word-boundary regex for a short acronym."""
    return re.compile(rf"\b{re.escape(term)}(?:s)?\b", re.IGNORECASE)


@dataclass(frozen=True)
class PainVocabulary:
    """Immutable pain-term vocabulary, grouped by kind.

    `context` holds anatomy and generic-biology terms (brain, cortex, CSF,
    inflammation). They are deliberately excluded from `all_terms` and
    `matches()`: a paper mentioning the brain is not thereby a pain study, and
    treating them as pain-defining would make auto-exclusion useless.
    """

    core: tuple[str, ...]
    conditions: tuple[str, ...]
    mechanisms: tuple[str, ...]
    context: tuple[str, ...] = ()

    @property
    def all_terms(self) -> tuple[str, ...]:
        """Pain-defining terms only, deduplicated, order preserved. Excludes context."""
        return tuple(dict.fromkeys(self.core + self.conditions + self.mechanisms))

    def query_string(self, joiner: str = " OR ", quote: bool = True) -> str:
        """Render the pain-defining vocabulary as a boolean query clause."""
        terms = [f'"{t}"' if quote else t for t in self.all_terms]
        return joiner.join(terms)

    def matches(self, text: str) -> tuple[str, ...]:
        """Return the pain-defining terms occurring in `text`, case-insensitively.

        Terms of length <= 4 (e.g. "CCI", "DRG") are matched as whole words
        (with an optional trailing "s") to avoid false positives like "CCI"
        inside "occipital". Longer terms keep the original substring match.
        """
        lowered = text.lower()
        found = []
        for t in self.all_terms:
            if len(t) <= SHORT_TERM_MAX_LEN:
                if _short_term_pattern(t).search(text):
                    found.append(t)
            elif t.lower() in lowered:
                found.append(t)
        return tuple(found)


def load_pain_vocabulary(path: Path | None = None) -> PainVocabulary:
    """Load the pain vocabulary from YAML.

    Args:
        path: Override for the config location. Defaults to
            conf/search/pain_vocabulary.yaml.

    Returns:
        A frozen PainVocabulary.

    Raises:
        FileNotFoundError: If the config file is absent.
        ValueError: If the vocabulary is empty.
    """
    vocab_path = path or DEFAULT_VOCAB_PATH
    if not vocab_path.exists():
        raise FileNotFoundError(f"Pain vocabulary not found: {vocab_path}")

    with open(vocab_path) as f:
        raw = yaml.safe_load(f) or {}

    vocab = PainVocabulary(
        core=tuple(raw.get("core", []) or []),
        conditions=tuple(raw.get("conditions", []) or []),
        mechanisms=tuple(raw.get("mechanisms", []) or []),
        context=tuple(raw.get("context", []) or []),
    )
    if not vocab.all_terms:
        raise ValueError(f"Pain vocabulary at {vocab_path} is empty")

    logger.info(
        "Loaded %d pain-defining terms (+%d context terms)",
        len(vocab.all_terms),
        len(vocab.context),
    )
    return vocab
