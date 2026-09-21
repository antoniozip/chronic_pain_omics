#!/usr/bin/env python
"""Re-resolve corpus-stated ChEBI assignments through the ontology, as a check.

`build_metabolite_id_map.py` trusts a depositing study when it supplies both a
name and a ChEBI accession (the `metabolights_maf` tier), because the submitter
knows what they measured. That trust is worth auditing: this re-resolves those
names independently through OLS4 and reports how often the two agree.

The result is a manuscript number, so it is written to a file the claims
registry reads rather than computed once by hand. It had drifted exactly the
way the registry exists to prevent — the manuscript said "The 255
corpus-derived assignments were re-resolved", which was true of the 502-row
map and became false when full-recall resolution took the corpus tier to 724.

Disagreements are expected and are not errors. They are ChEBI parent/child or
acid/conjugate-base pairs — thiamine against thiamine(1+), gluconolactone
against D-glucono-1,5-lactone — where the study and the ontology name the same
measurement at different levels of specificity.

Usage:
    python scripts/crosscheck_corpus_chebi.py --apply
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import logging
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MAP = REPO_ROOT / "conf" / "analysis" / "metabolite_id_map.csv"
CACHE = REPO_ROOT / "data" / "interim" / "metabolomics" / "ols4_chebi_cache.json"
OUT = REPO_ROOT / "results" / "meta" / "metabolomics" / "corpus_chebi_crosscheck.csv"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def _builder():
    """Load the map builder, to share its exact lookup and cache semantics."""
    path = REPO_ROOT / "scripts" / "build_metabolite_id_map.py"
    spec = importlib.util.spec_from_file_location("_map_builder", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write the result")
    args = ap.parse_args()

    mod = _builder()
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    with open(MAP) as f:
        rows = [r for r in csv.DictReader(f)
                if (r.get("evidence") or "").startswith("metabolights_maf")
                and r.get("chebi_id")]
    logger.info("%d corpus-stated assignments to re-resolve", len(rows))

    out, resolved, agree = [], 0, 0
    for i, r in enumerate(rows, 1):
        if i % 100 == 0:
            CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
            logger.info("  %d/%d (%d resolved, %d agree)", i, len(rows), resolved, agree)
        hit = mod.ols4_lookup(r["name"], cache)
        time.sleep(0.05)
        ols_id = hit[0] if hit else ""
        if ols_id:
            resolved += 1
            if ols_id == r["chebi_id"]:
                agree += 1
        out.append({"name": r["name"], "corpus_chebi_id": r["chebi_id"],
                    "ols4_chebi_id": ols_id,
                    "agrees": "" if not ols_id else str(ols_id == r["chebi_id"])})

    CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    pct = 100.0 * agree / resolved if resolved else 0.0
    logger.info("\n%d corpus-stated | %d also resolve through OLS4 | %d agree (%.1f%%)",
                len(rows), resolved, agree, pct)

    if not args.apply:
        logger.info("\nDry run. %d rows would be written to %s",
                    len(out), OUT.relative_to(REPO_ROOT))
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "corpus_chebi_id",
                                          "ols4_chebi_id", "agrees"])
        w.writeheader()
        w.writerows(out)
    logger.info("-> %s  %d rows", OUT.relative_to(REPO_ROOT), len(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
