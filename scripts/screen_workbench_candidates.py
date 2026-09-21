#!/usr/bin/env python
"""Screen the Metabolomics Workbench enumeration into a PRISMA candidate record.

Discovery for Workbench enumerates the whole repository (4,507 studies) rather
than querying it, because the keyword endpoint is not usable for systematic
retrieval — `pain` returns one study while `chronic fatigue` returns fifteen.
See `src/cp_multiomics/ingest/metabolomics_workbench.py`.

Screening therefore happens here, against one explicit expression, and this
script writes the record: every candidate that passed the automatic screen,
with the fields `literature/prisma/metabolomics_manual_review.csv` already
uses, so the two repositories' records are directly comparable.

**Verdicts are left to the reviewer.** A title carrying pain vocabulary is a
candidate, not an inclusion: psoriatic-arthritis microbiome profiling and
Schwann-cell disease modelling are in this list because of their words, not
their designs. `reason_code_proposed` is filled only where a documented title
pattern makes the exclusion unambiguous, exactly as the MetaboLights record
uses that column; `verdict` stays `pending` until a human decides.

The 4,485 studies excluded by the automatic screen are recorded as a count
rather than 4,485 identical rows: the screen is one regex over a re-fetchable
enumeration, so the exclusion is reproducible rather than asserted.

Usage:
    python scripts/screen_workbench_candidates.py            # dry run
    python scripts/screen_workbench_candidates.py --apply
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import logging
import re
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from cp_multiomics.ingest import metabolomics_workbench as wb  # noqa: E402

OUT = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_candidates.csv"
CACHE = REPO_ROOT / "data" / "raw" / "metabolomics" / "workbench_enumeration.json"
SUMMARY = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_counts.csv"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

FIELDS = ["accession", "species", "tier", "verdict", "reason_code",
          "reason_source", "reason_code_proposed", "platform", "year",
          "title", "url", "n_samples", "source"]

# Title patterns whose exclusion is unambiguous. Anything not matched here is
# left for the reviewer rather than guessed at — the same discipline that put
# `not_recorded` in the MetaboLights record instead of an invented code.
PROPOSED = (
    # A nerve disorder is not automatically a pain phenotype. Optic neuropathy
    # is a visual deficit; these studies carry no pain measure.
    (re.compile(r"(?i)optic neuropath"), "no_pain_phenotype"),
    # Cell-derived disease models have no patients and no case-control arm.
    # Anchored at the start of the stem only: a trailing \b after a truncated
    # stem never matches, so `\bhpsc\b` excludes "hPSCs" — the same defect that
    # cut a first screen of this repository from 23 candidates to 7.
    (re.compile(r"(?i)\b(hpsc|ipsc|cell line|organoid|in vitro)"),
     "no_pain_versus_control_contrast"),
    # Response-to-therapy designs compare responders with non-responders, or a
    # timepoint with itself. Neither is disease against unaffected control.
    (re.compile(r"(?i)(methotrexate (response|non-response)|following initiation"
                r"|normalization in|cognitive behavio(u)?ral therapy)"),
     "intervention_pre_post_not_case_control"),
)


def enumerate_studies(refresh: bool) -> list[dict]:
    """Fetch the repository, caching it so a re-screen costs nothing."""
    if CACHE.exists() and not refresh:
        studies = json.loads(CACHE.read_text())
        logger.info("using cached enumeration: %d studies (%s)",
                    len(studies), CACHE.relative_to(REPO_ROOT))
        return studies

    def _get(url: str, params: dict | None = None, timeout: int = 300):
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return json.load(r)

    studies = wb.fetch_all(_get)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(studies))
    return studies


def propose(title: str) -> str:
    for pattern, code in PROPOSED:
        if pattern.search(title):
            return code
    return ""


def rows_for(studies: list[dict], modality: str) -> tuple[list[dict], collections.Counter]:
    kept, excluded = wb.screen(studies, modality)
    rows = []
    for s in sorted(kept, key=lambda x: x["study_id"]):
        rec = wb.to_record(s, modality)
        rows.append({
            "accession": rec.accession,
            "species": rec.species[0] if rec.species else "",
            "tier": "prioritised",
            "verdict": "pending",
            "reason_code": "",
            "reason_source": "automatic_screen",
            "reason_code_proposed": propose(rec.title),
            "platform": rec.platform,
            "year": rec.year,
            "title": rec.title,
            "url": rec.url,
            "n_samples": rec.n_samples,
            "source": rec.source,
        })
    return rows, collections.Counter(e["reason_code"] for e in excluded)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write the records")
    ap.add_argument("--refresh", action="store_true", help="re-fetch the enumeration")
    args = ap.parse_args()

    studies = enumerate_studies(args.refresh)

    rows: list[dict] = []
    counts: list[dict] = []
    for modality in ("metabolomics", "lipidomics"):
        mod_rows, excluded = rows_for(studies, modality)
        rows.extend(mod_rows)
        logger.info("\n%s: %d enumerated, %d candidates, %d excluded",
                    modality, len(studies), len(mod_rows), sum(excluded.values()))
        for code, n in sorted(excluded.items()):
            logger.info("   %-24s %d", code, n)
        counts.append({
            "modality": modality,
            "source": "METABOLOMICS_WORKBENCH",
            "n_enumerated": len(studies),
            "n_candidates": len(mod_rows),
            **{f"excluded_{c}": n for c, n in sorted(excluded.items())},
        })
        for r in mod_rows:
            logger.info("   %-9s %-18s n=%5s %s %s", r["accession"],
                        r["species"][:18], r["n_samples"], r["title"][:58],
                        f"[{r['reason_code_proposed']}]" if r["reason_code_proposed"] else "")

    if not args.apply:
        logger.info("\nDry run. %d candidate rows would be written to %s",
                    len(rows), OUT.relative_to(REPO_ROOT))
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    count_fields = sorted({k for c in counts for k in c})
    with open(SUMMARY, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=count_fields)
        w.writeheader()
        w.writerows(counts)
    logger.info("\n-> %s  %d rows", OUT.relative_to(REPO_ROOT), len(rows))
    logger.info("-> %s", SUMMARY.relative_to(REPO_ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
