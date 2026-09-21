#!/usr/bin/env python3
"""Curated case/control assignment for studies the automatic screen cannot do.

`pipeline/auto_annotate_geo.py` handles a study whose contrast is one
characteristics column with a control-like value in it. The studies here are
the ones it cannot, and each is a different reason: the contrast is on the
sample title, or one arm has to be dropped, or the study deposits two designs
in one series. They are the bulk equivalent of `CONTRASTS` in
`pipeline/18_assign_single_cell_groups.py`, and reading these rules is how the
assignment gets reviewed, because there is no other record of it.

Every rule carries a note saying why it is not the obvious one. Samples a rule
does not name are dropped, and the dropped count is reported, so an arm cannot
shrink silently.

    python scripts/assign_rescued_groups.py --dry-run
    python scripts/assign_rescued_groups.py --apply
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE = REPO_ROOT / "data" / "raw" / "geo_cache"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


#: field: the characteristics key to read, or "title" for the sample title.
#: case / control: values assigned to each arm, lowercased and compared whole.
#: Anything else is dropped.
CURATED: dict[str, dict] = {
    "GSE269047": {
        "field": "diagnosis",
        "case": ["fm", "co-diagnosis"],
        "control": ["healthy control"],
        "note": ("Four arms: ME/CFS 8, fibromyalgia 10, co-diagnosed 16, "
                 "healthy 9. Fibromyalgia is in scope and myalgic "
                 "encephalomyelitis is not a pain diagnosis, so the ME/CFS-only "
                 "samples are dropped rather than entered as cases. The "
                 "co-diagnosed arm has fibromyalgia by definition."),
    },
    "GSE51981": {
        "field": "endometriosis severity",
        "case": ["moderate/severe", "minimal/mild"],
        "control": ["no uterine pelvic pathology"],
        "note": ("The automatic screen splits this on `endometriosis/no "
                 "endometriosis`, 77 against 71, and that control arm is not "
                 "pathology-free: 37 of its 71 carry `uterine pelvic "
                 "pathology`, which is fibroids or adenomyosis and is itself a "
                 "cause of pelvic pain. Over half the controls would be a pain "
                 "population, the same defect as GSE182983 in partial form. "
                 "The control arm is therefore restricted to the 34 samples "
                 "with no uterine pelvic pathology, giving 77 against 34. "
                 "Merging the two endometriosis severities into one case arm "
                 "is what the contrast audit flags and is not the problem: "
                 "both are endometriosis. Cycle phase stays unbalanced across "
                 "the arms and nothing downstream carries a covariate to "
                 "adjust for it."),
    },
    "GSE17504": {
        "field": "disease status",
        "case": ["endometriosis"],
        "control": ["no endometriosis"],
        "require": ("treatment group", {"control"}),
        "note": ("A 2x2: disease status crossed with 8-br-cAMP treatment, five "
                 "samples per cell. The treatment is given to half of each arm, "
                 "so leaving it in makes a drug effect part of the disease "
                 "effect. Restricted to the untreated cells, 5 against 5. Same "
                 "defect the single-cell arm handles for GSE254360."),
    },
    "GSE134056": {
        "field": "group",
        "case": ["disease"],
        "control": ["control"],
        "note": "Two arms, named plainly; here only because the value is the "
                "generic word 'disease' rather than a condition.",
    },
    "GSE7305": {
        "field": "title",
        "case_match": r"disease",
        "control_match": r"normal",
        "note": "The contrast is on the sample title, "
                "'endometrium/ovary-disease N' against 'endometrium-normal N'; "
                "no characteristics column carries it.",
    },
    # GSE310557 was here and is not a bulk study. Its deposit is a Matrix
    # Market triplet with `library source: transcriptomic single cell`, so its
    # seventeen "samples" are libraries. It had a groups file written into the
    # shared bulk cache before anyone looked, which is what library_is_single_cell()
    # now refuses. It belongs to the single-cell arm.
    "GSE11783": {
        "field": "disease state",
        "case": ["ic, ulcer"],
        "control": ["no"],
        "note": ("16 samples from 11 patients: the five interstitial cystitis "
                 "patients each deposit an ulcer and a non-ulcer biopsy. Both "
                 "would enter as independent cases and inflate n by five, so "
                 "only the lesional site is taken. 5 against 6."),
    },
    "GSE28242": {
        "field": "disease state",
        "case": ["pbs without lesions", "pbs with lesions"],
        "control": ["normal"],
        "note": ("PBS is painful bladder syndrome, not phosphate-buffered "
                 "saline, and the automatic screen reads it as the latter. "
                 "Both patient arms are cases: lesion-free patients are "
                 "patients. 8 against 5."),
    },
    "GSE58178": {
        "field": "tissure",
        "case": ["endometriotic tissue"],
        "control": ["normal endometrial tissue"],
        "note": "The characteristics key is misspelt 'tissure' in the deposit, "
                "which is also why it is not skipped as a tissue column.",
    },
    "GSE315857": {
        "field": "title",
        "case_match": r"^eosis",
        "control_match": r"^control",
        "note": "'eosis' is a truncation of endometriosis in the sample "
                "titles; no characteristics column carries the contrast.",
    },
    "GSE621": {
        "field": "title",
        "case_match": r"^ic patient",
        "control_match": r"^normal control",
        "note": ("Two designs in one series. Twelve samples are interstitial "
                 "cystitis and normal explants, a clean 6 against 6; the other "
                 "sixteen are cell lines given antiproliferative factor or "
                 "mock, which is a stimulus and not a diagnosis. Admitting the "
                 "series whole would enter a treatment effect as a disease "
                 "effect in more than half of it."),
    },
}


def parse(accession: str) -> tuple[list[str], list[str], dict[str, list[str]]]:
    """GSM accessions, titles, and characteristics field -> per-sample values."""
    paths = (sorted(CACHE.glob(f"{accession}-*_series_matrix.txt.gz"))
             + sorted(CACHE.glob(f"{accession}_series_matrix.txt.gz")))
    gsms: list[str] = []
    titles: list[str] = []
    chars: dict[str, list[str]] = {}
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("!series_matrix_table_begin"):
                    break
                if not line.startswith("!Sample_"):
                    continue
                cells = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                key, values = cells[0], [v.strip() for v in cells[1:]]
                if key == "!Sample_geo_accession":
                    gsms += values
                elif key == "!Sample_title":
                    titles += values
                elif key == "!Sample_characteristics_ch1":
                    names = [v.split(":", 1)[0].strip().lower()
                             for v in values if ":" in v]
                    if not names:
                        continue
                    field = max(set(names), key=names.count)
                    chars.setdefault(field, []).extend(
                        v.split(":", 1)[1].strip() if ":" in v else ""
                        for v in values)
    return gsms, titles, chars


def library_is_single_cell(accession: str) -> bool:
    """True when the deposit says its libraries are single cell.

    The reserved-accession guard in `auto_annotate_geo` only knows studies
    already listed as single-cell, so it cannot catch a new one. GSE310557 is
    "Expression profiling by high throughput sequencing" with a Matrix Market
    deposit and `library source: transcriptomic single cell`, and it was
    written a groups file in the shared bulk cache before anyone looked.
    """
    paths = (sorted(CACHE.glob(f"{accession}-*_series_matrix.txt.gz"))
             + sorted(CACHE.glob(f"{accession}_series_matrix.txt.gz")))
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("!series_matrix_table_begin"):
                    break
                if line.startswith("!Sample_library_source") and \
                        "single cell" in line.lower():
                    return True
    return False


def assign(accession: str, rule: dict) -> list[tuple[str, str]]:
    """(gsm, group) for every sample the rule names.

    `require` is a pre-filter on a second field, for a study whose contrast is
    crossed with something else. Without it a rule can only read one field, and
    a treatment crossed with the disease would enter as part of the disease
    effect.
    """
    if library_is_single_cell(accession):
        raise ValueError(
            f"{accession}: the deposit says its libraries are single cell. A "
            f"groups file in the shared GEO cache is what the bulk arm reads, "
            f"so this belongs in the single-cell arm instead.")
    gsms, titles, chars = parse(accession)
    if rule["field"] == "title":
        values = titles
    else:
        values = chars.get(rule["field"], [])
    if len(values) != len(gsms):
        raise ValueError(
            f"{accession}: field {rule['field']!r} has {len(values)} values "
            f"for {len(gsms)} samples")

    keep = set(gsms)
    if "require" in rule:
        field, allowed = rule["require"]
        req = chars.get(field, [])
        if len(req) != len(gsms):
            raise ValueError(f"{accession}: require field {field!r} has "
                             f"{len(req)} values for {len(gsms)} samples")
        keep = {g for g, v in zip(gsms, req) if v.strip().lower() in allowed}

    out = []
    for gsm, value in zip(gsms, values):
        if gsm not in keep:
            continue
        low = value.strip().lower()
        if "case_match" in rule:
            if re.search(rule["case_match"], low):
                out.append((gsm, "case"))
            elif re.search(rule["control_match"], low):
                out.append((gsm, "control"))
        elif low in rule["case"]:
            out.append((gsm, "case"))
        elif low in rule["control"]:
            out.append((gsm, "control"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="write the groups files; otherwise report only")
    args = ap.parse_args()

    for accession, rule in CURATED.items():
        pairs = assign(accession, rule)
        n_case = sum(1 for _, g in pairs if g == "case")
        n_ctrl = len(pairs) - n_case
        n_total = len(parse(accession)[0])
        logger.info("%-11s %2d case, %2d control, %2d dropped of %2d",
                    accession, n_case, n_ctrl, n_total - len(pairs), n_total)
        if n_case < 2 or n_ctrl < 2:
            raise ValueError(f"{accession}: an arm is below two samples")
        if args.apply:
            out = CACHE / f"{accession}_groups.csv"
            with out.open("w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["sample_id", "group"])
                writer.writerows(pairs)
    logger.info("%s %d curated study/studies",
                "wrote" if args.apply else "would write", len(CURATED))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
