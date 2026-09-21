#!/usr/bin/env python
"""Screen the retrieved GEO transcriptomic records into a PRISMA candidate record.

The transcriptomic arm never had an explicit screening record. Retrieval ran,
and a study entered the analysis if `auto_annotate_geo.py` could derive a
case/control split from its series matrix -- eligibility was decided implicitly
by whether the machinery worked, not by a recorded judgement. That was
tolerable at 110 records with one search; it is not, now that a disclosed
search amendment has taken the corpus to 283 and the reason a study is absent
matters.

This writes the record the arm should always have had: one row per retrieved
series, its verdict, and a reason code from the vocabulary the metabolomics and
proteomics screens already use.

**Screening is on titles.** The harmonized records carry no abstract, so this
reproduces the limitation of the Workbench screen: a title carrying pain
vocabulary is a candidate, not an inclusion, and the file check downstream is
what decides. Verdicts that need a human are left `pending` rather than
guessed.

**Phenotype scope is recorded, not decided here.**
`conf/analysis/phenotype_scope.csv` holds which conditions are in scope, when
that was decided and why, across every arm. Endometriosis is in scope
(confirmed 2026-08-28), which is what keeps 101 of the 165 candidates this
screen surfaces. That file also records an unresolved inconsistency: the
rheumatoid and juvenile idiopathic arthritis metabolomic cohorts were excluded
by scope while knee osteoarthritis and endometriosis are included.

**Studies already contributing to the analysis are marked, not re-judged.**
They were screened by the earlier process and their results are in the
manuscript; re-deciding them here would silently change published numbers.

Usage:
    python scripts/screen_transcriptomics_candidates.py            # dry run
    python scripts/screen_transcriptomics_candidates.py --apply
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import logging
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
HARMONIZED = REPO_ROOT / "data" / "interim" / "transcriptomics" / "harmonized.jsonl"
BASELINE = REPO_ROOT / "temp" / "baseline_tx" / "datasets.jsonl"
PER_STUDY = REPO_ROOT / "results" / "per_study" / "transcriptomics"
OUT = REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

SPECIES_KEEP = {"homo sapiens", "mus musculus", "rattus norvegicus"}

# Pain vocabulary: mechanisms and named conditions, matching the amended
# retrieval query. Stems are anchored at their start only -- a trailing \b
# after a truncated stem never matches, which is how `\barthrit\b` silently
# excluded "arthritis" from an earlier screen.
PAIN = re.compile(
    r"(?i)\b(pain|nocicept|noxious|analges|allodyn|hyperalges|neuralgi|migrain|"
    r"fibromyalg|arthrit|neuropath|sciatic|nerve injur|nerve ligat|neuroma|"
    r"complete freund|cfa\b|endometrios|cystitis|irritable bowel|dysmenorrh|"
    r"headache|crps|causalgi|chemotherapy.induced|cipn|low back|disc degener|"
    r"trigemin|temporomandib|postherpetic|spinal cord injur|dorsal root|drg\b)")

# Designs that cannot yield a patient-versus-control or model-versus-sham
# contrast, checked against the title.
NO_CONTRAST = re.compile(
    r"(?i)(ipsc|hpsc|organoid|cell line|in vitro|crispr|knock[- ]?in screen|"
    r"reprogramm|differentiation of|derived .{0,20}cells)")


def in_analysis() -> set[str]:
    """Accessions with a per-study effect table, i.e. already in the pool.

    Transcriptomic effect tables are flat files in one directory, unlike the
    metabolomic ones which sit in a directory per study. A glob written for the
    metabolomic layout matches nothing here and reports every pooled study as
    unscreened, which would mark 58 studies already in the manuscript as
    pending.
    """
    out = set()
    for path in glob.glob(str(PER_STUDY / "*_effects*.csv")):
        stem = os.path.basename(path)
        acc = stem.split("_effects")[0]
        if acc.startswith("GSE"):
            # Derived per-tissue ids (GSE241361_DRG) belong to their parent.
            out.add(acc.split("_")[0])
    return out


def screen(rec: dict) -> tuple[str, str]:
    """(verdict, reason_code) for one retrieved record."""
    title = rec.get("title") or ""
    species = rec.get("species_canonical") or []
    if isinstance(species, str):
        species = [species]
    if not any(s.strip().lower() in SPECIES_KEEP for s in species):
        return "exclude", "wrong_species"
    if not PAIN.search(title):
        return "exclude", "no_pain_phenotype"
    if NO_CONTRAST.search(title):
        return "exclude", "no_pain_versus_control_contrast"
    return "pending", "usable_pending_file_check"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write the record")
    args = ap.parse_args()

    recs = [json.loads(line) for line in HARMONIZED.open()]
    previously = ({json.loads(line)["accession"] for line in BASELINE.open()}
                  if BASELINE.exists() else set())
    pooled = in_analysis()
    logger.info("%d retrieved | %d retrieved before the amendment | %d already pooled",
                len(recs), len(previously), len(pooled))

    rows = []
    for rec in sorted(recs, key=lambda r: r["accession"]):
        acc = rec["accession"]
        species = rec.get("species_canonical") or []
        if isinstance(species, str):
            species = [species]
        if acc in pooled:
            verdict, reason = "include", "already_in_pooled_analysis"
        else:
            verdict, reason = screen(rec)
        rows.append({
            "accession": acc,
            "species": "; ".join(species),
            "n_samples": rec.get("n_samples", ""),
            "year": rec.get("year", ""),
            "retrieval": "original" if acc in previously else "amendment_2026_08_28",
            "verdict": verdict,
            "reason_code": reason,
            "reason_source": "automatic_screen",
            "title": rec.get("title", ""),
            "url": rec.get("url", ""),
        })

    by = collections.Counter((r["retrieval"], r["verdict"], r["reason_code"]) for r in rows)
    logger.info("\n%-22s %-9s %-38s %s", "retrieval", "verdict", "reason", "n")
    for (ret, verdict, reason), n in sorted(by.items()):
        logger.info("%-22s %-9s %-38s %d", ret, verdict, reason, n)

    new_pending = [r for r in rows
                   if r["retrieval"].startswith("amendment") and r["verdict"] == "pending"]
    logger.info("\n%d new candidates need a file check", len(new_pending))
    for r in sorted(new_pending, key=lambda r: -int(r["n_samples"] or 0))[:15]:
        logger.info("   %-11s n=%-5s %s", r["accession"], r["n_samples"], r["title"][:66])

    if not args.apply:
        logger.info("\nDry run. %d rows would be written to %s",
                    len(rows), OUT.relative_to(REPO_ROOT))
        return 0
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    logger.info("\n-> %s  %d rows", OUT.relative_to(REPO_ROOT), len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
