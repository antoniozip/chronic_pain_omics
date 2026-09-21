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

# Ordered: first match wins. Each entry is (canonical key, pattern).
_COHORT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("UKB", re.compile(r"uk ?biobank|\bukb\b", re.IGNORECASE)),
    ("FinnGen", re.compile(r"finngen", re.IGNORECASE)),
    ("23andMe", re.compile(r"23andme", re.IGNORECASE)),
]


def parse_sample_size(text: str) -> int:
    """Sum every comma-formatted integer in an initialSampleSize string.

    Returns -1 when the string contains no integer (unparseable).
    """
    matches = _INT_RE.findall(text or "")
    if not matches:
        return -1
    return sum(int(m.replace(",", "")) for m in matches)


def detect_cohort(text: str) -> str:
    """Return a canonical cohort key detected in the text, or "" if none."""
    for key, pattern in _COHORT_PATTERNS:
        if pattern.search(text or ""):
            return key
    return ""


def is_burden_study(trait: str) -> bool:
    """True if the trait names a gene-based burden analysis (not SNP-level GWAS)."""
    return bool(_BURDEN_RE.search(trait or ""))
