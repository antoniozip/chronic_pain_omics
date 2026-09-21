#!/usr/bin/env python
"""Pull the evidence needed to accept or reject each Workbench candidate.

`screen_workbench_candidates.py` screens 4,507 studies down to candidates on
title vocabulary alone, because the summary endpoint carries no abstract and no
design. A title is not an inclusion criterion, and screening on one gets both
kinds of error: ST002974 reads as a pain study and is RAW 264.7 cell pellets,
while ST003177 reads as a treatment-response study and carries 420 rheumatoid
arthritis samples against 438 healthy controls.

This script fetches, per candidate, the three things that actually decide it:

1. **`/factors`** — every sample's factor values. This is decisive and
   objective: a study with no unaffected level cannot yield a case-versus-
   control effect size at all, whatever its title says. It is the criterion
   that removed MTBLS9662 from the MetaboLights arm.
2. **`/metabolites`** — how many metabolites are named, and how many carry a
   `refmet_name`. RefMet is a standardised space that maps to ChEBI far more
   cleanly than free-text names, so this predicts how well a study would join
   the existing pool.
3. **`ST:STUDY_SUMMARY`** from `/mwtab/txt` — the abstract, for the phenotype
   judgement a human has to make.

**Verdicts are not written here.** The control-arm test is mechanical and its
result is reported as `no_control_arm` evidence; whether a rheumatoid arthritis
cohort or a diabetic neuropathy model counts as chronic pain for this study is
a scoping decision for the maintainer. This script assembles the evidence and
proposes; `verdict` stays `pending`.

Usage:
    python scripts/validate_workbench_candidates.py            # dry run
    python scripts/validate_workbench_candidates.py --apply
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import logging
import re
import sys
import time
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_candidates.csv"
OUT = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_validation.csv"
CACHE = REPO_ROOT / "data" / "interim" / "metabolomics" / "workbench_validation_cache.json"

REST = "https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/{ep}"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

SESSION = requests.Session()
TIMEOUT = 60

# mwtab carries the whole data table after the metadata header. Only the header
# is wanted, so the response is streamed and abandoned once the abstract is
# past — some of these studies are tens of megabytes.
MWTAB_HEADER_BYTES = 60_000

# Level values naming an unaffected arm. Anchored at the start of the stem, not
# wrapped in \b...\b: a trailing \b after a truncated stem never matches, which
# is how `\barthrit\b` silently excluded "arthritis" from a first screen.
#
# `hlt` and the general `non-<something>` clause were added after both produced
# false negatives on real studies: ST000218 names its healthy arm `HLT`, and
# ST001412 names its matched control `Obese non neuropathy` (44 v 44).
CONTROL_PATTERN = re.compile(
    r"(?i)\b(control|ctrl|healthy|health\b|\bhlt\b|normal|non[- ]?\w+|"
    r"sham|na(i|ï)ve|baseline|unaffected|wild[- ]?type|\bwt\b|vehicle|\bhc\b|"
    r"untreated|reference|lean)"
)

# Keys naming a timepoint rather than a group. A control-shaped level under one
# of these is a pre-treatment measurement of the *same* subjects, not an
# unaffected arm: ST001940 is `Timepoint = Baseline(34) / Post_CBT(34)`, 34 IBS
# patients measured twice, and reading `Baseline` as a control would pool a
# within-subject intervention as a case-versus-control effect.
TIMEPOINT_KEYS = re.compile(r"(?i)^(time ?point|time|visit|week|day|session|stage)s?$")

# Levels the `non-<something>` clause must NOT claim. A non-responder is a
# treated patient whose disease did not remit — the opposite of an unaffected
# control, and pooling one as a control would invert the contrast. ST004518 is
# `Response(45) / Non-response(15)`, which is a treatment-response design and
# carries no control arm at all.
NOT_CONTROL_PATTERN = re.compile(r"(?i)non[- ]?respon")

# Evidence that the samples are not organisms. Checked against factor values,
# where cell-line studies describe themselves even when the title does not.
#
# `in vitro` excludes `in vitro fertilisation`, which describes how a patient
# was treated, not what was assayed: ST000585 samples follicular fluid from
# women undergoing IVF and is a human cohort, not a cell culture.
INVITRO_PATTERN = re.compile(
    r"(?i)(cell pellet|cell line|raw 264|hek ?293|hela|sh[- ]?sy5y|"
    r"ipsc|hpsc|organoid|in vitro(?![ -]?fertili[sz])|cultured|passage \d)"
)


# Smallest arm that can yield an effect size, matching MIN_PER_GROUP in
# `pipeline/05j_metabolomics_da.py`. A "control level" of one sample is a
# label, not an arm.
MIN_ARM = 3


def is_design_factor(levels: collections.Counter, n_samples: int) -> bool:
    """Whether a factor key's levels are reused rather than per-sample.

    Used only to decide what to show a reviewer in `design_levels`; the
    control-arm test does not depend on it. Two earlier rules failed here and
    both are worth remembering:

    A hardcoded list of "specimen" key names was wrong in both directions —
    ST000585 puts its whole contrast in `Sample Type` (`Control(44)` against
    `Endometriosis patient(36)`) while ST002974 uses that key for 66 free-text
    cell-pellet descriptions. Replacing it with this reuse test then rejected
    ST003954, whose `Group1` groups the controls (`HC`, 10) but labels every
    case individually (`IBS-C-1` … `IBS-C-61`) — a genuine 10-against-61
    design wearing 62 levels.
    """
    n_levels = len(levels)
    if n_levels < 2:
        return False
    return n_levels <= max(2, n_samples / 3)


class FetchError(Exception):
    """A request that did not return data, as distinct from returning none."""


def get(acc: str, ep: str, cache: dict, attempts: int = 3) -> object:
    """Fetch one endpoint, memoised.

    **Failures are never cached.** A connection error cached as `None` is
    indistinguishable downstream from a study that genuinely has no
    metabolites, and would silently disqualify it — ST003177 (2,492 samples)
    reported 0 metabolites on the first run for exactly this reason. A failure
    raises instead, so the caller records it as missing evidence rather than as
    absent data.
    """
    key = f"{acc}/{ep}"
    if key in cache:
        return cache[key]
    last = None
    for attempt in range(attempts):
        try:
            r = SESSION.get(REST.format(acc=acc, ep=ep), timeout=TIMEOUT)
            r.raise_for_status()
            payload = r.json()
        except Exception as exc:
            last = exc
            time.sleep(1.0 * (attempt + 1))
            continue
        cache[key] = payload
        time.sleep(0.1)
        return payload
    logger.warning("  %s %s failed after %d attempts: %s", acc, ep, attempts, last)
    raise FetchError(f"{acc}/{ep}: {last}")


def rows_of(payload: object, id_field: str) -> list[dict]:
    """Normalise Workbench's two response shapes.

    A single-record response is the bare object, whose fields would otherwise
    be counted as records. Same trap as the summary endpoint.
    """
    if not payload:
        return []
    if isinstance(payload, dict):
        if id_field in payload:
            return [payload]
        return [v for v in payload.values() if isinstance(v, dict)]
    if isinstance(payload, list):
        return [v for v in payload if isinstance(v, dict)]
    return []


def parse_factors(rows: list[dict]) -> dict[str, collections.Counter]:
    """Split `Key:Value | Key:Value` strings into per-key level counts."""
    per_key: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for r in rows:
        for part in (r.get("factors") or "").split("|"):
            if ":" not in part:
                continue
            key, _, value = part.partition(":")
            key, value = key.strip(), value.strip()
            if key:
                per_key[key].update([value])
    return dict(per_key)


def has_control_arm(per_key: dict[str, collections.Counter],
                    n_samples: int) -> tuple[bool, str]:
    """Whether any design factor carries an unaffected level.

    Returns the matching level and its size, so a reviewer sees what the test
    fired on rather than trusting a bare boolean — which is how the `HLT` and
    `Obese non neuropathy` misses were caught.
    """
    for key, levels in per_key.items():
        if TIMEPOINT_KEYS.match(key.strip()):
            continue
        total = sum(levels.values())
        for level, n in levels.items():
            if NOT_CONTROL_PATTERN.search(level) or not CONTROL_PATTERN.search(level):
                continue
            # Both arms must be poolable. This is what separates a control arm
            # from a stray control-shaped label: ST002974's 66 cell-pellet
            # descriptions each cover one sample, so no level clears MIN_ARM,
            # while ST003954's oddly-encoded `HC(10)` against 61 individually
            # named cases clears it easily.
            if n >= MIN_ARM and (total - n) >= MIN_ARM:
                return True, f"{key}={level} (n={n} vs {total - n})"
    return False, ""


def looks_in_vitro(per_key: dict[str, collections.Counter], summary: str = "") -> str:
    """Evidence that the samples are not organisms.

    The abstract is searched as well as the factor values, because a study can
    describe its design nowhere else: ST003523 says "we used patient
    iPSC-derived motor neurons" while its factors read `Genotype = MUT / WT`,
    which is indistinguishable from an animal knockout.
    """
    for key, levels in per_key.items():
        for level in levels:
            if INVITRO_PATTERN.search(level):
                return f"{key}={level[:60]}"
    m = INVITRO_PATTERN.search(summary or "")
    if m:
        start = max(0, m.start() - 30)
        return "abstract: ..." + (summary[start:m.end() + 30]).strip()
    return ""


def study_summary(acc: str, cache: dict) -> str:
    key = f"{acc}/mwtab_summary"
    if key in cache:
        return cache[key] or ""
    text = ""
    try:
        with SESSION.get(REST.format(acc=acc, ep="mwtab/txt"),
                         timeout=TIMEOUT, stream=True) as r:
            r.raise_for_status()
            buf = r.raw.read(MWTAB_HEADER_BYTES, decode_content=True)
            text = buf.decode("utf-8", "replace")
    except Exception as exc:
        logger.warning("  %s mwtab failed: %s", acc, exc)
        return ""          # not cached: a failed fetch is not an empty abstract
    lines = [ln.split("\t", 1)[1].strip()
             for ln in text.splitlines()
             if ln.startswith("ST:STUDY_SUMMARY") and "\t" in ln]
    summary = " ".join(lines)
    cache[key] = summary
    time.sleep(0.1)
    return summary


def validate(acc: str, cache: dict) -> dict:
    errors = []
    try:
        fac_rows = rows_of(get(acc, "factors", cache), "local_sample_id")
    except FetchError as exc:
        fac_rows, _ = [], errors.append(f"factors: {exc}")
    per_key = parse_factors(fac_rows)
    control, control_evidence = has_control_arm(per_key, len(fac_rows))
    summary = study_summary(acc, cache)
    in_vitro = looks_in_vitro(per_key, summary)

    try:
        met_rows = rows_of(get(acc, "metabolites", cache), "metabolite_name")
    except FetchError as exc:
        met_rows, _ = [], errors.append(f"metabolites: {exc}")
    named = {r.get("metabolite_name", "") for r in met_rows if r.get("metabolite_name")}
    refmet = {r.get("refmet_name", "") for r in met_rows if r.get("refmet_name")}

    # A proposal must never rest on evidence that failed to arrive. Missing
    # evidence is reported as missing, not read as an absence.
    if errors:
        proposed = "evidence_incomplete"
    elif in_vitro:
        proposed = "no_pain_versus_control_contrast"
    elif not fac_rows:
        proposed = "factors_unavailable"
    elif not named:
        # No named metabolites and no MS_METABOLITE_DATA block in mwtab: the
        # Workbench counterpart of the five MetaboLights studies that declare
        # sample columns and leave every cell empty. Verified against mwtab for
        # ST000585 and ST003579, both of which lack the data block entirely
        # while known-good studies carry it.
        proposed = "no_abundances_deposited"
    elif not control:
        proposed = "no_control_arm_in_factors"
    else:
        proposed = ""

    return {
        "accession": acc,
        "fetch_errors": "; ".join(errors),
        "n_samples_with_factors": len(fac_rows),
        "factor_keys": "; ".join(sorted(per_key)),
        "design_levels": " | ".join(
            f"{k}={','.join(f'{lv}({n})' for lv, n in c.most_common(6))}"
            for k, c in per_key.items()
            if is_design_factor(c, len(fac_rows)))[:400],
        "has_control_arm": control,
        "control_evidence": control_evidence,
        "in_vitro_evidence": in_vitro,
        "n_metabolites": len(named),
        "n_refmet": len(refmet),
        "reason_code_proposed": proposed,
        "study_summary": summary[:1200],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write the validation table")
    ap.add_argument("--accessions", help="comma-separated subset")
    args = ap.parse_args()

    if args.accessions:
        accs = [a.strip() for a in args.accessions.split(",")]
    else:
        with open(CANDIDATES) as f:
            accs = [r["accession"] for r in csv.DictReader(f)]

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    results = []
    for i, acc in enumerate(accs, 1):
        logger.info("[%d/%d] %s", i, len(accs), acc)
        results.append(validate(acc, cache))
        CACHE.write_text(json.dumps(cache))

    logger.info("\n%-10s %6s %6s %7s  %-34s %s", "study", "n", "mets",
                "refmet", "control", "proposed")
    for r in results:
        logger.info("%-10s %6d %6d %7d  %-34s %s", r["accession"],
                    r["n_samples_with_factors"], r["n_metabolites"],
                    r["n_refmet"], (r["control_evidence"] or "NONE")[:34],
                    r["reason_code_proposed"])

    n_ok = sum(1 for r in results if not r["reason_code_proposed"])
    logger.info("\n%d of %d clear every mechanical gate "
                "(control arm, not in vitro, abundances deposited)", n_ok, len(results))

    if not args.apply:
        logger.info("\nDry run. %d rows would be written to %s",
                    len(results), OUT.relative_to(REPO_ROOT))
        return 0

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    logger.info("\n-> %s  %d rows", OUT.relative_to(REPO_ROOT), len(results))
    return 0


if __name__ == "__main__":
    sys.exit(main())
