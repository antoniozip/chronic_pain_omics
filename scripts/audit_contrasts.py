#!/usr/bin/env python3
"""Check every study's case/control split against what its GEO metadata says.

Why this exists
---------------
`groups.csv` declares which samples are cases and which are controls, and
nothing downstream ever questions it. A study that splits on drug treatment,
on microarray channel, on cell type or on timepoint produces a perfectly
well-formed effect table, enters the pool, and reports a treatment effect as a
pain effect.

That is not hypothetical. The 2026-08-29 audit of the search amendment's 46
studies disqualified six on exactly this: GSE248593 split Cy5 against Cy3
channels of the same six samples, GSE296883 entered six drug-treatment arms as
"case", GSE343056 made baclofen-*treated* diabetic animals the case arm and
dropped the untreated one, GSE318478 collapsed four diet arms into 13 v 7.
The same check had never been run over the 55 studies of the original
retrieval, which are the ones carrying the published CCI/SNI/SNL/CFA/CIPN
strata.

How it works
------------
For each study, every ``!Sample_characteristics_ch1`` field is parsed into
key -> value per sample and joined to ``groups.csv`` by GSM accession. A field
that **perfectly separates** the arms -- no value of it appearing in both --
is the field the contrast is really on. That is exactly right when the field
describes the disease or the model (``surgery type: Sham`` v ``SNI``), and is
a finding when it describes anything else.

Fields are therefore classified, and a study is reported when a separating
field is not a model field. Reported, not judged: whether a given split is
wrong is a reading of the study, which is why the output is a review sheet
rather than a verdict.

Two structural checks run alongside, for the other two failure shapes the
amendment audit found:

  * **channel duplication** -- a two-colour array deposits each sample twice,
    so n is half what the arm counts say (GSE248593);
  * **pseudoreplication** -- a subject/animal/donor field with fewer distinct
    values than samples means the same individual appears more than once
    (GSE72428: 19 samples from 9 baboons; GSE313775: 22 subjects x 3 T-cell
    subsets entered as 66).

Usage:
    python scripts/audit_contrasts.py
    python scripts/audit_contrasts.py --accessions GSE15041,GSE18803
    python scripts/audit_contrasts.py --out temp/contrast_audit.csv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
import sys
from collections import Counter, defaultdict
from functools import reduce
from math import gcd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CACHE = REPO_ROOT / "data" / "raw" / "geo_cache"
PRISMA = REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv"

logger = logging.getLogger(__name__)

#: Field names that legitimately carry the case/control contrast. A field here
#: separating the arms is the design working, not a finding.
MODEL_FIELD = re.compile(
    r"surgery|injury|model|group|condition|disease|status|diagnosis|phenotype"
    r"|treatment group|cci|sni|snl|cfa|sham|naive|pain|control|patient|case"
    r"|genotype/variation|cohort|subject group|experimental",
    re.I,
)

#: Field names that describe something other than the phenotype. A field here
#: separating the arms means the contrast is on that instead.
CONFOUND_FIELD = re.compile(
    r"channel|dye|cy3|cy5|label\b"                    # two-colour arrays
    r"|drug|dose|compound|vehicle|inhibitor|agonist|antagonist"
    r"|diet|feed|nutrition"
    r"|time|day\b|week\b|hour|timepoint|age"
    r"|cell type|celltype|cell subset|subset|sorted|fraction|population"
    r"|batch|scan|run\b|slide|array|chip"
    r"|tissue|region|site|organ"
    r"|sex|gender|strain",
    re.I,
)

TREATMENT_VALUE = re.compile(
    r"vehicle|saline|treated|treatment|mg/kg|µg|ug/ml|dose|drug|inhibitor"
    r"|agonist|antagonist|baclofen|morphine|gabapentin|paclitaxel|cisplatin",
    re.I,
)

#: Per-subject covariates. A person's BMI does not change between their own
#: aliquots, so if every BMI value appears exactly three times there are three
#: samples per subject. This is the one signal that catches both shapes of
#: inflated n at once: GSE313775's 22 subjects x 3 T-cell subsets (BMI repeats
#: 3x) and GSE248593's two-colour array (age repeats 2x), neither of which
#: separates the arms and so neither of which a confound check can see.
COVARIATE_FIELD = re.compile(
    r"\bage\b|bmi|body mass|weight|height|race|ethnic|education|smok"
    r"|disease_stage|disease stage|duration", re.I,
)

#: Channel tokens a two-colour array appends to otherwise identical titles.
CHANNEL_TOKEN = re.compile(
    r"[\s_-]*\b(cy\s?[357]|alexa\s?\d+|635|532|channel\s?[12ab]|ch[12])\b",
    re.I,
)

SUBJECT_FIELD = re.compile(
    r"subject|animal|patient|donor|individual|participant|mouse id|rat id"
    r"|\bid\b|identifier", re.I,
)


def series_matrices(acc: str) -> list[Path]:
    """Every series matrix cached for an accession, not just the first.

    A SuperSeries whose first matrix is a methylation half hides its expression
    SubSeries from any check that reads one file, which is how GSE102937 and
    GSE51296 survived the 2026-08-28 screen.
    """
    return sorted(CACHE.glob(f"{acc}-*_series_matrix.txt.gz")) + \
        sorted(CACHE.glob(f"{acc}_series_matrix.txt.gz"))


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace")


def parse_matrix(path: Path) -> tuple[list[str], dict[str, list[str]], list[str]]:
    """GSM accessions in order, characteristics field -> per-sample values, titles."""
    gsms: list[str] = []
    titles: list[str] = []
    chars: dict[str, list[str]] = defaultdict(list)
    raw_char_lines: list[list[str]] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            cells = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key, values = cells[0], cells[1:]
            if key == "!Sample_geo_accession":
                gsms = [v.strip() for v in values]
            elif key == "!Sample_title":
                titles = [v.strip() for v in values]
            elif key == "!Sample_characteristics_ch1":
                raw_char_lines.append([v.strip() for v in values])

    # Each characteristics line is one field. The field name is the part before
    # the first colon, which GEO repeats on every sample; take the most common
    # so a single malformed cell does not rename the field.
    for row in raw_char_lines:
        names = [c.split(":", 1)[0].strip().lower()
                 for c in row if ":" in c]
        if not names:
            continue
        field = max(set(names), key=names.count)
        vals = [c.split(":", 1)[1].strip() if ":" in c else "" for c in row]
        if field in chars:                       # same field twice: keep both
            field = f"{field} (2)"
        chars[field] = vals
    return gsms, dict(chars), titles


def load_groups(acc: str) -> dict[str, str]:
    path = CACHE / f"{acc}_groups.csv"
    if not path.exists():
        return {}
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("sample_id") or "").strip()
            grp = (row.get("group") or "").strip().lower()
            if sid and grp in ("case", "control"):
                out[sid] = grp
    return out


def separating_fields(chars: dict[str, list[str]], gsms: list[str],
                      groups: dict[str, str]) -> list[tuple[str, str]]:
    """Fields no value of which appears in both arms, with their level split."""
    out = []
    for field, values in chars.items():
        if len(values) != len(gsms):
            continue
        by_value: dict[str, set[str]] = defaultdict(set)
        for gsm, val in zip(gsms, values):
            arm = groups.get(gsm)
            if arm and val:
                by_value[val.lower()].add(arm)
        if len(by_value) < 2:
            continue
        if all(len(arms) == 1 for arms in by_value.values()):
            levels = "; ".join(
                f"{v}={next(iter(a))}" for v, a in sorted(by_value.items())[:6])
            out.append((field, levels))
    return out


def audit(acc: str) -> dict | None:
    groups = load_groups(acc)
    mats = series_matrices(acc)
    if not groups or not mats:
        return None

    gsms: list[str] = []
    chars: dict[str, list[str]] = {}
    titles: list[str] = []
    for path in mats:
        g, c, t = parse_matrix(path)
        if not gsms:
            gsms, chars, titles = g, dict(c), t
        elif g == gsms:                          # same samples, more fields
            for k, v in c.items():
                chars.setdefault(k, v)
        else:                                    # a second platform's samples
            gsms += g
            titles += t
            for k, v in c.items():
                chars[k] = chars.get(k, [""] * (len(gsms) - len(g))) + v

    # Most older deposits carry the contrast in the sample title and nothing
    # else -- 18 of the original 55 declare no characteristics field that
    # separates the arms at all, and two declare no characteristics at all.
    # Treating the title as a field is what lets the model-field test mean
    # anything on those studies.
    if titles and len(titles) == len(gsms):
        chars.setdefault("sample title", titles)

    n_case = sum(1 for g in groups.values() if g == "case")
    n_ctrl = len(groups) - n_case
    sep = separating_fields(chars, gsms, groups)

    flags: list[str] = []
    detail: list[str] = []
    for field, levels in sep:
        if CONFOUND_FIELD.search(field) and not MODEL_FIELD.search(field):
            flags.append("separates_on_" + re.sub(r"\W+", "_", field)[:40])
            detail.append(f"{field} -> {levels}")
        elif TREATMENT_VALUE.search(levels) and not MODEL_FIELD.search(field):
            flags.append("treatment_values_in_" + re.sub(r"\W+", "_", field)[:34])
            detail.append(f"{field} -> {levels}")

    # An arm drawing on several levels of a separating categorical field is a
    # merger of design groups. Sometimes that is right -- sham and naive are
    # both controls -- and sometimes it hides a second variable: GSE318478's
    # case arm merges `soy oil diet_injury` with `fish oil diet_cci`, so the
    # contrast carries diet as well as injury. Reported with the levels so the
    # reviewer can tell the two apart; this cannot be decided mechanically.
    for field, levels in sep:
        if field == "sample title" or COVARIATE_FIELD.search(field):
            continue
        per_arm: dict[str, set[str]] = defaultdict(set)
        for value, arm in (pair.split("=", 1)[::-1]
                           for pair in levels.split("; ") if "=" in pair):
            per_arm[value].add(arm)
        if any(len(v) >= 2 for v in per_arm.values()):
            flags.append("arms_merge_levels_of_"
                         + re.sub(r"\W+", "_", field)[:32])
            detail.append(f"{field} -> {levels}")

    # Duplicated sample titles, once the dye label is stripped: a two-colour
    # array deposits each biological sample twice, so n is half what the arm
    # counts say. GSE248593's titles differ only by a trailing "Cy5"/"Cy3".
    assigned = [CHANNEL_TOKEN.sub("", t).strip().lower()
                for gsm, t in zip(gsms, titles) if gsm in groups]
    assigned = [a for a in assigned if a]
    if assigned and len(set(assigned)) < len(assigned):
        flags.append("duplicate_sample_titles")
        detail.append(f"{len(assigned)} assigned samples, "
                      f"{len(set(assigned))} distinct titles once the dye "
                      f"label is stripped")

    # A per-subject covariate repeating in an exact multiple: every subject
    # contributes that many samples.
    for field, values in chars.items():
        if not COVARIATE_FIELD.search(field) or len(values) != len(gsms):
            continue
        vals = [v for gsm, v in zip(gsms, values) if gsm in groups and v]
        counts = Counter(vals)
        if len(counts) < 3:
            continue
        factor = reduce(gcd, counts.values())
        if factor >= 2:
            flags.append(f"n_inflated_{factor}x")
            detail.append(f"{field}: {len(counts)} distinct values, each "
                          f"appearing a multiple of {factor} times, so "
                          f"{len(vals)} samples are ~{len(vals) // factor} subjects")

    # A subject field with fewer levels than samples: repeated individuals.
    for field, values in chars.items():
        if not SUBJECT_FIELD.search(field) or len(values) != len(gsms):
            continue
        vals = [v for gsm, v in zip(gsms, values) if gsm in groups and v]
        if vals and len(set(vals)) * 1.5 < len(vals):
            flags.append("repeated_subjects_in_" + re.sub(r"\W+", "_", field)[:34])
            detail.append(f"{field}: {len(vals)} samples, {len(set(vals))} subjects")

    # Reported as a column, not a flag. A deposit whose characteristics do not
    # encode the contrast is under-annotated, which is common and is not
    # evidence that the contrast is wrong; flagging it drowned the seven
    # findings that matter under 47 that do not.
    model_sep = any(MODEL_FIELD.search(f) for f, _ in sep)

    return {
        "accession": acc, "n_case": n_case, "n_control": n_ctrl,
        "n_fields": len(chars), "n_separating": len(sep),
        "model_field_separates": "yes" if model_sep else "no",
        "separating_fields": " | ".join(f for f, _ in sep),
        "flags": ";".join(sorted(set(flags))),
        "detail": " || ".join(detail)[:600],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", help="comma-separated subset")
    ap.add_argument("--wave", choices=["original", "amendment", "all"],
                    default="all")
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "literature" / "prisma"
                    / "transcriptomics_contrast_audit.csv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    with open(PRISMA, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["verdict"] == "include"]
    if args.wave == "original":
        rows = [r for r in rows if r.get("retrieval") == "original"]
    elif args.wave == "amendment":
        rows = [r for r in rows if r.get("retrieval") != "original"]
    accs = [r["accession"] for r in rows]
    if args.accessions:
        want = set(args.accessions.split(","))
        accs = [a for a in accs if a in want]

    results, missing = [], []
    for acc in accs:
        res = audit(acc)
        if res is None:
            missing.append(acc)
            continue
        results.append(res)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]) if results else
                           ["accession"])
        w.writeheader()
        w.writerows(results)

    flagged = [r for r in results if r["flags"]]
    logger.info("audited %d studies, %d without cached metadata",
                len(results), len(missing))
    if missing:
        logger.info("  no groups.csv or series matrix: %s", ", ".join(missing))
    logger.info("%d flagged -> %s", len(flagged), args.out)
    for r in sorted(flagged, key=lambda x: x["accession"]):
        logger.info("  %-12s %s", r["accession"], r["flags"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
