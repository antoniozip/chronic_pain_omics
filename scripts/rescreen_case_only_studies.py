#!/usr/bin/env python3
"""Re-screen studies the automatic screen called case-only, from their metadata.

Why this exists
---------------
75 human studies carry `no_case_control_split_in_series_matrix`, and **every
one of those verdicts was set by `auto_annotate_geo` with no manual review**.
That screen looks for the vocabulary in `conf/analysis/default.yaml`, whose
`case_labels` are pain terms -- "chronic pain", "neuropathic", "fibromyalgia",
"CRPS" -- and whose `control_labels` are "healthy", "control", "normal".

A study of a disease the vocabulary does not name therefore reports no split
even when its samples are labelled as plainly as `disease: endometriosis`
against `disease: control`. The absence looks identical to a study that
genuinely has no control arm. GSE276193 is the precedent for what that costs:
the same class of automatic verdict, and the series matrix showed 8 cases
against 10 controls once anyone fetched it.

What it does
------------
For each accession it parses every cached series matrix and, for each
characteristics field, asks whether the field's values split the samples into
a control-like group and a non-control-like group. Control-likeness is lexical
and deliberately broad, because the point is to surface candidates a
pain-specific vocabulary cannot see.

**It reports, it does not decide.** A field that separates samples may be
describing cycle phase or tissue rather than disease, which is a reading of the
study. The output is a review sheet; changing a verdict is a separate act.

Usage:
    python scripts/rescreen_case_only_studies.py
    python scripts/rescreen_case_only_studies.py --accessions GSE11691,GSE23339
    python scripts/rescreen_case_only_studies.py --out temp/rescreen.csv
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import audit_contrasts as ac  # noqa: E402

CANDIDATES = REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv"
DEFAULT_OUT = REPO_ROOT / "literature" / "prisma" / "transcriptomics_rescreen_audit.csv"
TARGET_REASON = "no_case_control_split_in_series_matrix"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

#: Control-like values. Broader than the pain vocabulary in
#: conf/analysis/default.yaml, which is what missed these studies, but not
#: broader than a diagnosis: `without` and `no` negate a *diagnosis* here,
#: never a feature of one.
#: Bare `without` read "PBS without lesions" in GSE28242 as a control and put
#: five painful bladder syndrome patients in the control arm -- a well-formed
#: split that is wrong, which is the failure this whole audit exists to catch.
_NEGATED = r"(?:disease|diseased|endometrio\w*|pain\w*|neuropath\w*|cystitis|" \
           r"symptoms?|diagnosis|cancer|fibromyalgia)"
CONTROL_RE = re.compile(
    r"\b(?:control|healthy|normal|unaffected|negative|sham|na[iï]ve|ctrl|ctl|"
    r"vehicle|baseline|reference"
    r"|non[- ]?" + _NEGATED
    + r"|(?:no|without)\s+" + _NEGATED
    + r"|" + _NEGATED + r"[- ]free)\b", re.I)

#: Fields that separate samples without describing disease state. A field
#: matching this is reported with a flag rather than proposed as the contrast.
#: Short tokens are anchored. Unanchored `age` matches inside "stage", which
#: flagged GSE141549's "disease stage" as non-disease and let a cycle-phase
#: field outrank it on the arm's largest study, 408 samples.
NON_DISEASE_RE = re.compile(
    r"cycle|phase|menstrual|tissue|cell ?type|region|timepoint|time ?point|"
    r"passage|batch|platform|\bage\b|\bsex\b|gender|\bbmi\b|\brace\b|"
    r"ethnic|treatment|drug|\bdose\b|channel|replicate|donor|subject|"
    r"patient ?id|library|\brun\b", re.I)


def is_control(value: str) -> bool:
    return bool(value) and bool(CONTROL_RE.search(value))


def split_on(values: list[str]) -> tuple[int, int, str, str]:
    """Control count, non-control count, and a label for each side."""
    controls = [v for v in values if is_control(v)]
    others = [v for v in values if v and not is_control(v)]
    ctrl_label = Counter(v.lower() for v in controls).most_common(1)
    case_label = Counter(v.lower() for v in others).most_common(1)
    return (len(controls), len(others),
            case_label[0][0] if case_label else "",
            ctrl_label[0][0] if ctrl_label else "")


def candidate_fields(chars: dict[str, list[str]], n: int) -> list[dict]:
    """Every field that puts at least one sample on each side of control."""
    found = []
    for field, values in chars.items():
        if len(values) != n:
            continue
        n_ctrl, n_case, case_label, ctrl_label = split_on(values)
        if n_ctrl == 0 or n_case == 0:
            continue
        found.append({
            "field": field,
            "n_case": n_case,
            "n_control": n_ctrl,
            "case_label": case_label,
            "control_label": ctrl_label,
            "n_values": len({v.lower() for v in values if v}),
            "looks_non_disease": bool(NON_DISEASE_RE.search(field)),
        })
    # A disease-looking field with few levels is the likeliest contrast.
    found.sort(key=lambda f: (f["looks_non_disease"], f["n_values"]))
    return found


def rescreen(accession: str) -> dict:
    """One review row for one accession."""
    matrices = ac.series_matrices(accession)
    if not matrices:
        return {"accession": accession, "status": "no_series_matrix_cached"}

    gsms: list[str] = []
    chars: dict[str, list[str]] = {}
    titles: list[str] = []
    for path in matrices:
        part_gsms, part_chars, part_titles = ac.parse_matrix(path)
        offset = len(gsms)
        gsms += part_gsms
        titles += part_titles
        for field, values in part_chars.items():
            # Pad so fields present in only one matrix of a SuperSeries still
            # line up with the samples they belong to.
            chars.setdefault(field, [""] * offset)
            chars[field] += values
    for field in chars:
        chars[field] += [""] * (len(gsms) - len(chars[field]))
    if titles and len(titles) == len(gsms):
        chars["sample title"] = titles

    fields = candidate_fields(chars, len(gsms))
    row = {
        "accession": accession,
        "status": "candidate_contrast_found" if fields else "no_control_like_value",
        "n_samples": len(gsms),
        "n_fields": len(chars),
        "n_candidate_fields": len(fields),
    }
    if fields:
        best = fields[0]
        row.update({
            "field": best["field"],
            "n_case": best["n_case"],
            "n_control": best["n_control"],
            "case_label": best["case_label"],
            "control_label": best["control_label"],
            "field_looks_non_disease": "yes" if best["looks_non_disease"] else "no",
            "other_candidate_fields": "; ".join(
                f["field"] for f in fields[1:6]),
        })
    return row


COLUMNS = ["accession", "status", "n_samples", "n_fields", "n_candidate_fields",
           "field", "n_case", "n_control", "case_label", "control_label",
           "field_looks_non_disease", "other_candidate_fields"]


def targets() -> list[str]:
    with CANDIDATES.open() as fh:
        return [r["accession"] for r in csv.DictReader(fh)
                if "Homo" in (r.get("species") or "")
                and r.get("reason_code") == TARGET_REASON]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", help="comma-separated, default: every target")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    accessions = (args.accessions.split(",") if args.accessions else targets())
    rows = [rescreen(a) for a in accessions]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    tally = Counter(r["status"] for r in rows)
    for status, n in tally.most_common():
        logger.info("%5d  %s", n, status)
    disease = sum(1 for r in rows
                  if r.get("field_looks_non_disease") == "no")
    logger.info("%5d  of those, on a field that does not look non-disease", disease)
    logger.info("wrote %s (%d rows)", args.out, len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
