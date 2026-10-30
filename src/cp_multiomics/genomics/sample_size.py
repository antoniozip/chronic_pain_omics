"""Pure parsers for GWAS Catalog study metadata.

The GWAS Catalog stores sample size as free text (`initialSampleSize`) and has
no machine-readable cohort field, so both must be recovered by parsing. These
functions drive cohort-independence dedup (see plan/cohort_expansion_plan.md
§7.7), so they carry the heaviest test burden in this stage.
"""

from __future__ import annotations

import re

_INT_RE = re.compile(r"\b(\d{1,3}(?:,\d{3})*|\d+)\b")
_BURDEN_RE = re.compile(r"\bburden\b|gene-based", re.IGNORECASE)

# Numbers that are not sample sizes: a UKB field id, an age range, a cohort
# count. `initialSampleSize` free text carries all three ("UKB data field 3773",
# "aged 40-69", "from 3 cohorts"), and summing them inflated N — which is the
# METAL SAMPLESIZE weight and the cohort tie-break, so the error is a weight
# error, not a cosmetic one.
_NON_SAMPLE_RE = re.compile(
    r"(?:data\s+field|aged)\s+\d[\d,]*(?:\s*[-\u2013]\s*\d[\d,]*)?"
    r"|from\s+\d[\d,]*\s+cohorts?",
    re.IGNORECASE,
)

# A parenthetical breakdown of a total stated just before it.
_PAREN_RE = re.compile(r"\(([^()]*)\)")

# Ordered: first match wins. Each entry is (canonical key, pattern).
_COHORT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("UKB", re.compile(r"uk ?biobank|\bukb\b", re.IGNORECASE)),
    ("FinnGen", re.compile(r"finngen", re.IGNORECASE)),
    ("23andMe", re.compile(r"23andme", re.IGNORECASE)),
]


def _sample_size_tokens(text: str) -> list[tuple[int, int]]:
    """(start, value) for every integer that could be a sample size."""
    ignored = [m.span() for m in _NON_SAMPLE_RE.finditer(text)]
    out: list[tuple[int, int]] = []
    for m in _INT_RE.finditer(text):
        if any(a <= m.start() < b for a, b in ignored):
            continue
        out.append((m.start(), int(m.group(1).replace(",", ""))))
    return out


def parse_sample_size(text: str) -> int:
    """Return the sample size stated in an initialSampleSize string.

    Prefers a standalone total over the parenthetical breakdown that repeats
    it ("12,345 ... (6,789 cases and 5,556 controls)" is 12,345, not 24,690),
    and ignores integers introduced by "data field", "aged" or "from N
    cohorts". Returns -1 when no sample size can be recovered.
    """
    text = text or ""
    tokens = _sample_size_tokens(text)
    if not tokens:
        return -1

    for paren in _PAREN_RE.finditer(text):
        inner = _sample_size_tokens(paren.group(1))
        if not inner:
            continue
        before = [value for start, value in tokens if start < paren.start()]
        if before and sum(value for _, value in inner) == before[-1]:
            return before[-1]

    return sum(value for _, value in tokens)


def detect_cohort(text: str) -> str:
    """Return a canonical cohort key detected in the text, or "" if none."""
    for key, pattern in _COHORT_PATTERNS:
        if pattern.search(text or ""):
            return key
    return ""


def is_burden_study(trait: str) -> bool:
    """True if the trait names a gene-based burden analysis (not SNP-level GWAS)."""
    return bool(_BURDEN_RE.search(trait or ""))
