#!/usr/bin/env python
"""Resolve MetaboLights metabolite names to ChEBI, into conf/analysis/metabolite_id_map.csv.

The usable metabolomics studies deposit their metabolites in three identifier
namespaces — ChEBI, HMDB, and Metabolon/Workbench internal ids — so nothing is
measured in three studies in any single namespace and `meta.min_studies = 3`
yields an empty correction family. Every study does carry a *name*, so names are
the only space in which the arm can meet. This maps that space onto ChEBI.

It is the metabolomic counterpart of `05i_peptide_to_gene.py`, and the mapping
lives in a reviewable CSV for the same reason: assigning a name to a compound is
a chemical judgement, not a constant.

**Evidence has three tiers, and the strongest is the corpus itself.**

1. ``metabolights_maf`` — a study deposits both a name and a ChEBI id in its own
   MAF. That is the submitter's own assignment for this exact measurement, and
   it needs no external lookup. It also resolves names ChEBI does not carry as a
   label or synonym at all: "alpha-ketoisovaleric acid" fails every OLS4 query
   but MTBLS2774 states CHEBI:16530 for it.
2. ``ols4_label`` — the name equals a ChEBI term's primary label, compared
   case-insensitively.
3. ``ols4_synonym`` — the name equals one of a ChEBI term's exact synonyms.
   This is what recovers "D-Malic acid" as CHEBI:30796, whose label is
   "(R)-malic acid".

OLS4's ``exact=true`` is not exact: querying "L-tyrosine" returns
"L-tyrosinate(1-)" and "methyl L-tyrosinate" among its hits. The equality test
is applied here rather than trusted to the service.

**Stereochemistry is preserved.** Names are normalised by lower-casing and
dropping non-alphanumerics, which keeps the ``l-``/``d-`` prefix as part of the
token, so L-malic and D-malic stay separate compounds. An earlier estimate of
the reachable set stripped those prefixes and reported 170 metabolites at
k>=3; preserving them gives 162, and the difference is 8 pairs of enantiomers
that are not the same molecule.

**Resolve every name, not only the shared ones.** This script originally
resolved names measured in at least two studies, on the reasoning that a name
appearing once cannot contribute to a pooled estimate whatever it maps to.
That reasoning is wrong, and this script's own output is the counterexample:
seven ChEBI ids here carry two names each — Erucamide and 13Z-Docosenamide are
both CHEBI:142245, Azelaic acid and Nonanedioic acid are both CHEBI:48131 —
and CHEBI:48131 reaches the pooled k=3 family. A join is created by the
*accessions* agreeing, not by the spellings agreeing. Under the old filter
those seven surfaced only because both spellings happened to be shared
already; the identical case where each spelling appears in one study apiece
was never looked up, and 1,982 of the 1,988 names still keyed ``NAME:`` in the
per-study tables are exactly that. ``--shared-only`` restores the old
behaviour, so the 502-row map remains reproducible.

Usage:
    python scripts/build_metabolite_id_map.py            # dry run, all names
    python scripts/build_metabolite_id_map.py --apply
    python scripts/build_metabolite_id_map.py --shared-only   # pre-2026-08-27
"""

from __future__ import annotations

import argparse
import collections
import csv
import io
import json
import logging
import re
import sys
import time
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
REVIEW = REPO_ROOT / "literature" / "prisma" / "metabolomics_manual_review.csv"
WB_REVIEW = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_candidates.csv"
MWTAB = "https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/mwtab/txt"
OUT = REPO_ROOT / "conf" / "analysis" / "metabolite_id_map.csv"
CACHE = REPO_ROOT / "data" / "interim" / "metabolomics" / "ols4_chebi_cache.json"

FILES_API = "https://www.ebi.ac.uk/metabolights/ws/studies/{acc}/files"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{acc}/{name}"
OLS4 = "https://www.ebi.ac.uk/ols4/api/search"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

MIN_STUDIES = 2   # below this a name cannot pool, whatever it resolves to

UNNAMED = re.compile(
    r"^(cluster[_\- ]?\d+|unknown.*|unidentified.*|peak[_\- ]?\d+|m\d+t\d+"
    r"|feature[_\- ]?\d+|metabolite[_\- ]?\d+|compound[_\- ]?\d+|\d+(\.\d+)?)$",
    re.I,
)


def norm(name: str) -> str:
    """Lower-case and strip non-alphanumerics, keeping l-/d- stereo prefixes."""
    return re.sub(r"[^a-z0-9]+", "", name.strip().lower())


# One connection, reused. A fresh TLS handshake per lookup is most of the cost
# of ~5,700 of them.
SESSION = requests.Session()

# 120 s was the original timeout, chosen when only 548 shared names were
# resolved and a slow one cost little. Resolving every name walks into
# pathological inputs (see `name_variants`), and at 120 s apiece a few hundred
# of those dominate the whole run. A query OLS4 can answer returns in well
# under a second; one it cannot will not become answerable by waiting.
TIMEOUT = 15

# Flush the cache this often. Without it the cache is written only after the
# whole loop, so a run has no progress signal and is not resumable — a killed
# run at 2h17m discarded every lookup it had made.
FLUSH_EVERY = 200

# Queries longer than this are not sent. OLS4 does not answer a full IUPAC
# string with bracket notation — it times out rather than returning no match,
# so each costs the whole timeout and cannot succeed anyway. Roughly 50 such
# names accounted for most of a 2 h 17 m run.
#
# The threshold is set from the data, not guessed: across the 496 assignments
# in the pre-2026-08-27 map, the longest name **ever resolved through OLS4**
# is 43 characters and the 99th percentile is 35. (The 59-character maximum
# over all tiers is a corpus-stated accession, which needs no lookup.) 72
# leaves generous headroom for long lipid nomenclature such as
# "PC(17:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))" while still skipping the ~99-character
# offenders. An earlier cap of 100 would have skipped none of them — the string
# that triggered this is 99 characters.
#
# A skipped name is recorded as unresolved, never dropped.
MAX_QUERY_LEN = 72

# When true, `_ols4_query` answers only from the cache and never touches the
# network. OLS4 rate-limits a sustained client hard — measured at roughly one
# lookup per 45 s after a few thousand requests, against ~1 s from a cold
# start — so a full pass can stall for hours with thousands of names already
# resolved and stranded, because the map is written only after the loop. This
# lets the run be finished from what is already known: names never attempted
# are recorded as unresolved, exactly as a failed lookup is, and a later run
# picks them up once the limit resets.
CACHE_ONLY = False


def get_json(url: str, **params) -> dict:
    r = SESSION.get(url, params=params or None, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def name_variants(name: str) -> list[str]:
    """Candidate queries for one MAF name, best first.

    MAF rows do not always hold a compound name. Some hold a pipe-separated
    synonym blob, all of it in the `metabolite_identification` field:

        "5-methoxy-2,2-dimethyl-...-8-one|5-Methoxy-...|...|GUT-70"

    Queried whole, that matches nothing and times out; its last field is the
    compound name that would have resolved. Each part is therefore tried in
    turn. Parts are ordered shortest first because a trivial name resolves
    where a full IUPAC string does not, and the exact label/synonym equality
    test downstream is what keeps that from being a loose match.

    The old `>= min_studies` filter never met these: a vendor string that long
    does not collide across studies, so it was skipped rather than handled.
    """
    raw = name.strip()
    parts = [p.strip() for p in raw.split("|")] if "|" in raw else [raw]
    # dict.fromkeys keeps first-seen order among equal-length parts, so the
    # result is deterministic and the map reproducible.
    ordered = sorted(dict.fromkeys(p for p in parts if p), key=len)
    return [p for p in ordered if len(p) <= MAX_QUERY_LEN]


def collect_maf_names(accessions: list[str]) -> tuple[dict, dict]:
    """Return (normalised name -> studies) and (normalised name -> {chebi: studies})."""
    studies: dict[str, set[str]] = collections.defaultdict(set)
    corpus: dict[str, dict[str, set[str]]] = collections.defaultdict(
        lambda: collections.defaultdict(set))
    display: dict[str, str] = {}
    for acc in accessions:
        files = get_json(FILES_API.format(acc=acc)).get("study", []) or []
        for fname in sorted(f["file"] for f in files
                            if f["file"].startswith("m_") and f["file"].endswith(".tsv")):
            text = requests.get(FTP.format(acc=acc, name=fname), timeout=120).content.decode(
                "utf-8", "replace")
            rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
            if not rows:
                continue
            header = [c.strip() for c in rows[0]]
            di = header.index("database_identifier") if "database_identifier" in header else -1
            mi = (header.index("metabolite_identification")
                  if "metabolite_identification" in header else -1)
            if mi < 0:
                continue
            for r in rows[1:]:
                if len(r) <= mi:
                    continue
                raw = r[mi].strip()
                if not raw or UNNAMED.match(raw):
                    continue
                key = norm(raw)
                if not key:
                    continue
                studies[key].add(acc)
                display.setdefault(key, raw)
                if di >= 0 and len(r) > di and r[di].strip().upper().startswith("CHEBI"):
                    corpus[key][r[di].strip().upper()].add(acc)
    return {"studies": studies, "display": display}, corpus


def collect_workbench_names(accessions: list[str]) -> tuple[dict, dict]:
    """Metabolite names from the included Metabolomics Workbench studies.

    Returns the same (studies, display) shape as `collect_maf_names`, so both
    repositories feed one name corpus and `n_studies` counts across them. A
    name measured in a MetaboLights study and a Workbench study is one name
    seen twice, which is the entire point of resolving to a shared accession.

    **RefMet names are added alongside the study's own.** RefMet is a
    standardised naming space, so it resolves to ChEBI where a vendor
    abbreviation cannot: ST003984 deposits `Ac-Orn`, which no ontology carries,
    beside the RefMet name `N-Acetylornithine`, which resolves cleanly. Both
    are recorded because `05j` keys on whichever the study emits, and it
    prefers RefMet where present.
    """
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from cp_multiomics.metabolomics.workbench_mwtab import parse_mwtab

    studies: dict[str, set[str]] = collections.defaultdict(set)
    display: dict[str, str] = {}
    for acc in accessions:
        try:
            text = requests.get(MWTAB.format(acc=acc), timeout=300).content.decode(
                "utf-8", "replace")
        except Exception as exc:
            logger.warning("  %s mwtab failed: %s", acc, exc)
            continue
        names: set[str] = set()
        for analysis in parse_mwtab(text):
            names.update(analysis.values)
            names.update(analysis.refmet.values())
        for raw in names:
            if not raw or UNNAMED.match(raw):
                continue
            key = norm(raw)
            if not key:
                continue
            studies[key].add(acc)
            display.setdefault(key, raw)
        logger.info("  %s: %d names", acc, len(names))
    return studies, display


def ols4_lookup(name: str, cache: dict) -> tuple[str, str, str] | None:
    """Return (chebi_id, chebi_label, tier) for an exact label or synonym match.

    Cached under the MAF name as given, so a pipe blob is looked up once
    however many variants it expands to.
    """
    if name in cache:
        hit = cache[name]
        return tuple(hit) if hit else None
    for variant in name_variants(name):
        hit = _ols4_query(variant, cache)
        if hit:
            cache[name] = list(hit)
            return hit
    cache[name] = None
    return None


def _ols4_query(name: str, cache: dict) -> tuple[str, str, str] | None:
    """One OLS4 exact-match query, memoised under the queried string.

    That string is a *variant* when the caller expanded a pipe blob, so a
    variant shared by two blobs costs one request, not two. `ols4_lookup` then
    records the result under the original MAF name as well.

    OLS4's `exact=true` is not exact — querying "L-tyrosine" returns
    "L-tyrosinate(1-)" and "methyl L-tyrosinate" — so equality is tested here
    rather than trusted to the service. That test is also what makes trying
    several variants safe: a variant matches only if a ChEBI label or exact
    synonym equals it outright.
    """
    if CACHE_ONLY:
        return None
    try:
        data = get_json(OLS4, q=name, ontology="chebi", rows=10,
                        fieldList="obo_id,label,synonym", queryFields="label,synonym")
    except Exception as exc:
        logger.warning("  OLS4 failed for %r: %s", name[:80], exc)
        return None
    docs = data.get("response", {}).get("docs", []) or []
    target = name.strip().lower()
    for doc in docs:                      # primary label wins over synonym
        if (doc.get("label") or "").strip().lower() == target and doc.get("obo_id"):
            cache[name] = [doc["obo_id"], doc["label"], "ols4_label"]
            return doc["obo_id"], doc["label"], "ols4_label"
    for doc in docs:
        for syn in doc.get("synonym") or []:
            if syn.strip().lower() == target and doc.get("obo_id"):
                cache[name] = [doc["obo_id"], doc["label"], "ols4_synonym"]
                return doc["obo_id"], doc["label"], "ols4_synonym"
    cache[name] = None
    return None


def chebi_label(chebi_id: str, cache: dict) -> str:
    """Look up a ChEBI term's primary label by id, for review legibility."""
    key = "id:" + chebi_id
    if key in cache:
        return cache[key] or ""
    if CACHE_ONLY:
        return ""
    try:
        data = get_json(OLS4, q=chebi_id, ontology="chebi", rows=5,
                        fieldList="obo_id,label", queryFields="obo_id")
    except Exception:
        cache[key] = ""
        return ""
    for doc in data.get("response", {}).get("docs", []) or []:
        if doc.get("obo_id") == chebi_id:
            cache[key] = doc.get("label", "")
            return cache[key]
    cache[key] = ""
    return ""


def select_names(studies: dict[str, set[str]], min_studies: int,
                 shared_only: bool) -> dict[str, set[str]]:
    """Choose which names to resolve.

    Defaults to every name. A name seen in a single study is still worth an
    accession, because the join that matters is between accessions: two
    studies spelling one compound differently meet only after both spellings
    resolve, and neither spelling is shared. ``shared_only`` restores the
    pre-2026-08-27 filter for reproducing the 502-row map.
    """
    if not shared_only:
        return dict(studies)
    return {k: v for k, v in studies.items() if len(v) >= min_studies}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write the map")
    ap.add_argument("--shared-only", action="store_true",
                    help="resolve only names measured in >= --min-studies "
                         "studies (the pre-2026-08-27 behaviour)")
    ap.add_argument("--min-studies", type=int, default=MIN_STUDIES,
                    help="threshold for --shared-only; ignored otherwise")
    ap.add_argument("--from-cache-only", action="store_true",
                    help="resolve only from the OLS4 cache, making no network "
                         "lookups; names never attempted are recorded as "
                         "unresolved. Use to finish a run that OLS4 rate-limiting "
                         "has stalled, then re-run normally once it resets.")
    args = ap.parse_args()

    global CACHE_ONLY
    CACHE_ONLY = args.from_cache_only
    if CACHE_ONLY:
        logger.info("cache-only: no OLS4 requests will be made\n")

    with open(REVIEW) as f:
        accs = [r["accession"] for r in csv.DictReader(f) if r["verdict"] == "include"]
    logger.info("Reading MAFs from %d included MetaboLights studies", len(accs))
    names, corpus = collect_maf_names(accs)
    studies, display = names["studies"], names["display"]

    if WB_REVIEW.exists():
        with open(WB_REVIEW) as f:
            wb_accs = [r["accession"] for r in csv.DictReader(f)
                       if r["verdict"] == "include"]
        logger.info("Reading mwtab from %d included Workbench studies", len(wb_accs))
        wb_studies, wb_display = collect_workbench_names(wb_accs)
        for key, accs_seen in wb_studies.items():
            studies[key].update(accs_seen)
            display.setdefault(key, wb_display[key])

    n_shared = sum(1 for v in studies.values() if len(v) >= args.min_studies)
    shared = select_names(studies, args.min_studies, args.shared_only)
    logger.info("%d distinct names, %d measured in >=%d studies; resolving %d\n",
                len(studies), n_shared, args.min_studies, len(shared))

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    rows, unresolved, ambiguous = [], [], []
    order = sorted(shared, key=lambda k: (-len(studies[k]), k))
    started = time.time()
    for done, key in enumerate(order, 1):
        # Flush periodically: the cache is the run's only durable progress, and
        # writing it once at the end makes a long run both unobservable and
        # unresumable.
        if done % FLUSH_EVERY == 0:
            CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
            rate = done / max(time.time() - started, 1e-9)
            logger.info("  %d/%d names (%.1f/s, %d resolved, ~%.0f min left)",
                        done, len(order), rate, len(rows),
                        (len(order) - done) / rate / 60)
        name = display[key]
        chebi = label = tier = ""
        # Tier 1: the corpus states it.
        if key in corpus:
            ids = corpus[key]
            if len(ids) > 1:
                # Two studies gave one name two ChEBI ids. In this corpus these
                # are ChEBI parent/child or acid/conjugate-base pairs
                # (thiamine vs thiamine(1+), gluconolactone vs
                # D-glucono-1,5-lactone), so picking one is a chemical
                # judgement about which entity the measurement refers to.
                # Guessing would either merge two entities or split one, so the
                # row is kept with an empty chebi_id and the candidates
                # recorded, exactly as peptide_gene_map.csv keeps its
                # unassigned targets on the record.
                ambiguous.append((name, sorted(ids)))
                rows.append({
                    "name": name,
                    "normalised_name": key,
                    "chebi_id": "",
                    "chebi_label": "",
                    "n_studies": len(studies[key]),
                    "studies": ";".join(sorted(studies[key])),
                    "evidence": "ambiguous_in_corpus:" + "|".join(sorted(ids)),
                })
                continue
            chebi = next(iter(ids))
            tier = "metabolights_maf:" + ",".join(sorted(next(iter(ids.values()))))
        else:
            hit = ols4_lookup(name, cache)
            time.sleep(0.05)
            if hit:
                chebi, label, tier = hit
        if not chebi:
            unresolved.append((name, len(studies[key])))
            continue
        rows.append({
            "name": name,
            "normalised_name": key,
            "chebi_id": chebi,
            "chebi_label": label,
            "n_studies": len(studies[key]),
            "studies": ";".join(sorted(studies[key])),
            "evidence": tier,
        })

    # Backfill labels for corpus-tier rows. A reviewer cannot audit a bare
    # CHEBI:16919, and the corpus supplies an id without a name for it.
    labels = {r["chebi_id"]: r["chebi_label"] for r in rows if r["chebi_label"]}
    for r in rows:
        if r["chebi_id"] and not r["chebi_label"]:
            if r["chebi_id"] not in labels:
                labels[r["chebi_id"]] = chebi_label(r["chebi_id"], cache)
                time.sleep(0.05)
            r["chebi_label"] = labels[r["chebi_id"]]

    CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))

    assigned = [r for r in rows if r["chebi_id"]]
    by_tier = collections.Counter(r["evidence"].split(":")[0] for r in assigned)
    logger.info("resolved %d of %d names", len(assigned), len(shared))
    for t, c in sorted(by_tier.items()):
        logger.info("   %-20s %d", t, c)
    logger.info("unresolved: %d | ambiguous within corpus: %d", len(unresolved), len(ambiguous))

    # What actually matters: ChEBI ids now measured in >= 3 studies.
    per_chebi: dict[str, set[str]] = collections.defaultdict(set)
    for r in assigned:
        per_chebi[r["chebi_id"]].update(r["studies"].split(";"))
    dist = collections.Counter(len(v) for v in per_chebi.values())
    logger.info("\nChEBI ids by study count: %s", {k: dist[k] for k in sorted(dist)})
    logger.info("   in >=3 studies: %d", sum(v for k, v in dist.items() if k >= 3))

    if ambiguous:
        logger.info("\nambiguous (one name, several ChEBI ids in the corpus):")
        for name, ids in ambiguous[:10]:
            logger.info("   %-40s %s", name[:40], ids)
    if unresolved:
        logger.info("\ntop unresolved names by study count:")
        for name, k in sorted(unresolved, key=lambda x: -x[1])[:12]:
            logger.info("   k=%d  %s", k, name[:64])

    if not args.apply:
        logger.info("\nDry run. %d rows would be written to %s",
                    len(rows), OUT.relative_to(REPO_ROOT))
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["name", "normalised_name", "chebi_id", "chebi_label",
              "n_studies", "studies", "evidence"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    logger.info("\n-> %s  %d rows", OUT.relative_to(REPO_ROOT), len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
