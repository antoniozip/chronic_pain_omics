"""Assign each GSM of the single-cell arm to a case or control arm.

Why this exists
---------------
`literature/prisma/single_cell_contrast_audit.csv` counts keyword hits over
`!Sample_characteristics` and `!Sample_title`. Those are *field* hits, so the
audit can say a contrast is present and cannot say which sample belongs to
which arm. Everything downstream needs the second thing, and it cannot be had
mechanically: four of the twelve studies here carry the contrast in a field no
keyword list would find, and five carry samples that must be excluded for
reasons that have nothing to do with the contrast.

So the rules are curated, per study, in `CONTRASTS` below, and each one
records why. Reading them is how the arm's group assignment gets reviewed;
there is no other record of it.

What the rules can express
--------------------------
Each arm is one (field, regex) pair. `field` names a `!Sample_characteristics`
key, or `title` for `!Sample_title`, or `chars` for every characteristics value
of that sample joined -- which is what GSE198608 needs, its two submission
batches having spelled the same field two different ways.

`require` narrows the study to the samples eligible for the contrast at all,
`subject` names the biological unit so that two libraries of one animal are
summed rather than counted twice, and `drop` records a sample removed for a
stated reason. Every sample the rules do not place in an arm is written to the
PRISMA record with the reason it was not placed, so nothing leaves the arm
silently.

Outputs
-------
    data/interim/single_cell/groups/{accession}_groups.csv
        sample_id, group, subject, title -- read by step 19 and step 05.
    literature/prisma/single_cell_group_assignments.csv
        every GSM of all 18 series, placed or not, with the deciding reason.
    literature/prisma/single_cell_candidates.csv
        the screening record step 05 reads, keyed by accession.

Usage:
    python pipeline/18_assign_single_cell_groups.py
    python pipeline/18_assign_single_cell_groups.py --accession GSE134003
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE = REPO_ROOT / "data" / "raw" / "geo_cache"
REGISTRY = REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv"
PRISMA = REPO_ROOT / "literature" / "prisma"
GROUPS_DIR = REPO_ROOT / "data" / "interim" / "single_cell" / "groups"

#: A study needs this many samples in each arm. Mirrors `da.min_samples_per_group`
#: in conf/analysis/default.yaml, which is what step 05 enforces; applying it
#: here means a study that cannot be analysed is refused with a reason rather
#: than failing silently two steps later.
MIN_PER_GROUP = 2

#: Per-study contrast rules. See the module docstring for the grammar.
CONTRASTS: dict[str, dict] = {
    "GSE134003": {
        "case": ("title", r"^SNI\d+$"),
        "control": ("title", r"^SI\d+$"),
        "note": (
            "SI is superficial injury, the study's sham surgery control, and is "
            "in no control vocabulary -- which is why the keyword audit reported "
            "zero control-like fields for a study whose design statement reads "
            "'Transcriptomes from 10 SI surgery mice and 10 SNI mice'. One SNI "
            "sample was withdrawn by the submitters for DRG contamination, so "
            "the series carries 19."
        ),
    },
    "GSE155622": {
        "require": [("technology", r"^10x genomics$")],
        "case": ("treatment", r"^SNI \S+$"),
        "control": ("treatment", r"^untreated$"),
        "note": (
            "296 of the 314 GSMs are one Smart-seq2 *neuron* each, not one "
            "animal each, and the metadata carries no animal id, so they cannot "
            "be aggregated to a biological replicate at all. Treating them as "
            "samples is the cell-as-replicate error this arm exists to avoid, so "
            "the study enters on its 18 10x libraries only. Those span both "
            "platform records, which is why the filter is on technology rather "
            "than on the series matrix file."
        ),
    },
    "GSE162807": {
        "require": [("title", r"^Mouse ")],
        "subject": ("title", r"^Mouse microglia (\d+)"),
        "case": ("surgery", r"^SNI$"),
        "control": ("surgery", r"^(Sham|Naive)$"),
        "note": (
            "The three human samples are all naive, so the human half carries no "
            "contrast and the study is mouse-only -- which also keeps it out of "
            "the multi-species exclusion that would otherwise drop it from the "
            "species-resolved pools. Three animals were sequenced twice (103 and "
            "103-Redo, 106 and 106-Redo, 307A and 307B); the libraries are summed "
            "into one sample per animal rather than entered twice."
        ),
    },
    "GSE179640": {
        "require": [("method", r"^scRNA-seq$"),
                    ("sample location", r"^Eutopic$"),
                    ("tissue", r"^Endometrium$")],
        "case": ("condition", r"^Endometriosis$"),
        "control": ("condition", r"^Control$"),
        "note": (
            "The series mixes 24 bulk RNA-seq libraries into the same record, "
            "and its ectopic peritoneum and ovary samples have no control "
            "counterpart -- a lesion has nothing to be compared against in a "
            "woman without lesions. The contrast the bulk arm already runs on "
            "this disease is eutopic endometrium, patient against control, so "
            "that is the contrast taken here. Organoid and cell-hashing "
            "libraries are excluded by the same filter."
        ),
    },
    "GSE198608": {
        "case": ("chars", r"SNC \d"),
        "control": ("chars", r"Sham"),
        "note": (
            "Two submission batches spell the same field 'post-surgical time "
            "point' and 'post-surgical time_point', and the earlier batch omits "
            "a cell-type field so every later field is shifted by one. Matching "
            "the joined characteristics rather than a named key survives both. "
            "The seven GSM78392xx samples are all sham."
        ),
    },
    "GSE214411": {
        "case": ("disease", r"^endometriosis$"),
        "control": ("disease", r"^Control$"),
        "note": (
            "The three GPL11154 samples carry no disease field and are dropped "
            "for that reason, leaving the ten GPL24676 samples that do. Cycle "
            "phase is balanced across the arms."
        ),
    },
    "GSE253345": {
        "case": ("title", r"^CCI_.*_GEX$"),
        "control": ("title", r"^Control.*_GEX$"),
        "note": (
            "Half the series is ATAC, which is chromatin accessibility and not "
            "expression, and one sample is Smart-seq3. That leaves four gene "
            "expression libraries: one CCI and three control. One case is below "
            "da.min_samples_per_group, so the study is recorded as excluded "
            "rather than analysed -- and the exclusion is on sample count, not "
            "on the contrast, which is real."
        ),
    },
    "GSE254360": {
        "require": [("stimulation", r"^(Light touch|no stim)$")],
        "case": ("injury state", r"^SNI$"),
        "control": ("injury state", r"^Uninjured$"),
        "note": (
            "Stimulation is crossed with injury but not balanced: three "
            "uninjured animals received a hot water drop and no SNI animal did. "
            "Keeping them would put a stimulus with no counterpart entirely in "
            "the control arm, so they are dropped and the contrast runs 3v3 "
            "within each of the two stimuli that are matched."
        ),
    },
    "GSE283177": {
        "case": ("treatment", r"^PSNL surgery$"),
        "control": ("treatment", r"^sham surgery$"),
        "note": "Partial sciatic nerve ligation against sham, 3v3, balanced across two batches.",
    },
    "GSE289659": {
        "case": ("sample group", r"^Injured$"),
        "control": ("sample group", r"^(Control|Sham)$"),
        "note": (
            "Sham and naive are both controls, as elsewhere in this corpus. The "
            "field labelled genotype is sex; all four injured animals are female "
            "and the controls are 8 female and 4 male, so sex does not separate "
            "the arms."
        ),
    },
    "GSE328175": {
        "case": ("treatment", r"^SNI$"),
        "control": ("treatment", r"^Sham$"),
        "note": (
            "The two SNI+Transplant animals are a treatment arm and are dropped. "
            "Entering treated animals as cases is the error the 2026-08-29 "
            "contrast audit disqualified GSE343056 for."
        ),
    },
}

#: Studies the contrast audit found no usable case/control split in, and the
#: one whose deposit turned out not to be single-cell at all. Listed so the
#: screening record covers all 18 rather than only the ones that pass.
NO_CONTRAST = {
    "GSE203191": "no_case_control_split_in_series_matrix",
    "GSE213216": "all_samples_case_like_no_control_arm",
    "GSE222182": "all_samples_case_like_no_control_arm",
    "GSE266026": "all_samples_case_like_no_control_arm",
    "GSE266265": "all_samples_case_like_no_control_arm",
    "GSE276193": "bulk_rna_seq_not_single_cell_in_this_deposit",
    "GSE312510": "all_samples_case_like_no_control_arm",
}


def series_matrices(accession: str) -> list[Path]:
    """Every cached series matrix for an accession, platform splits included."""
    return sorted(CACHE.glob(f"{accession}-*_series_matrix.txt.gz")) + \
        sorted(CACHE.glob(f"{accession}_series_matrix.txt.gz"))


def parse_matrix(path: Path) -> list[dict]:
    """One record per GSM: accession, title, and each characteristics field."""
    gsms: list[str] = []
    titles: list[str] = []
    char_lines: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            cells = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key, values = cells[0], [v.strip() for v in cells[1:]]
            if key == "!Sample_geo_accession":
                gsms = values
            elif key == "!Sample_title":
                titles = values
            elif key == "!Sample_characteristics_ch1":
                char_lines.append(values)

    records = [{"sample_id": g, "title": titles[i] if i < len(titles) else "",
                "_chars": []} for i, g in enumerate(gsms)]
    for row in char_lines:
        if len(row) != len(gsms):
            continue
        # GEO repeats the field name on every sample; take the most common so a
        # single malformed cell cannot rename the field.
        names = [c.split(":", 1)[0].strip().lower() for c in row if ":" in c]
        if not names:
            continue
        field = max(set(names), key=names.count)
        for i, cell in enumerate(row):
            value = cell.split(":", 1)[1].strip() if ":" in cell else cell
            records[i].setdefault(field, value)
            records[i]["_chars"].append(value)
    for rec in records:
        rec["chars"] = " | ".join(rec.pop("_chars"))
    return records


def field_value(record: dict, field: str) -> str:
    return str(record.get(field, ""))


def matches(record: dict, spec: tuple[str, str]) -> bool:
    field, pattern = spec
    return bool(re.search(pattern, field_value(record, field), re.I))


def assign(accession: str, rule: dict) -> tuple[list[dict], str]:
    """Place every GSM of one study, returning the rows and a status."""
    matrices = series_matrices(accession)
    if not matrices:
        return [], "no_series_matrix_cached"

    records: list[dict] = []
    seen: set[str] = set()
    for path in matrices:
        for rec in parse_matrix(path):
            if rec["sample_id"] in seen:
                continue
            seen.add(rec["sample_id"])
            records.append(rec)

    subject_spec = rule.get("subject")
    rows: list[dict] = []
    for rec in records:
        row = {"accession": accession, "sample_id": rec["sample_id"],
               "title": rec["title"], "group": "", "subject": rec["sample_id"],
               "reason": ""}

        failed = [f for f, p in rule.get("require", []) if not matches(rec, (f, p))]
        if failed:
            row["reason"] = "excluded_not_eligible_on_" + \
                "_and_".join(re.sub(r"\W+", "_", f) for f in failed)
        elif matches(rec, rule["case"]):
            row["group"] = "case"
            row["reason"] = "case_on_" + re.sub(r"\W+", "_", rule["case"][0])
        elif matches(rec, rule["control"]):
            row["group"] = "control"
            row["reason"] = "control_on_" + re.sub(r"\W+", "_", rule["control"][0])
        else:
            row["reason"] = "excluded_matches_neither_arm"

        if subject_spec and row["group"]:
            found = re.search(subject_spec[1], field_value(rec, subject_spec[0]), re.I)
            if found:
                row["subject"] = f"{accession}_{found.group(1)}"
        rows.append(row)

    return rows, "assigned"


def study_verdict(rows: list[dict]) -> tuple[str, str, int, int]:
    """Screening verdict for one study from its assigned rows."""
    subjects: dict[str, str] = {}
    for row in rows:
        if row["group"]:
            subjects[row["subject"]] = row["group"]
    counts = Counter(subjects.values())
    n_case, n_control = counts.get("case", 0), counts.get("control", 0)
    if n_case < MIN_PER_GROUP or n_control < MIN_PER_GROUP:
        return ("exclude",
                f"below_min_samples_per_group_{n_case}_case_{n_control}_control",
                n_case, n_control)
    return "include", "usable_case_control_contrast", n_case, n_control


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", action="append",
                        help="restrict to these accessions (repeatable)")
    args = parser.parse_args()

    registry = pd.read_csv(REGISTRY).set_index("accession")
    wanted = args.accession or list(registry.index)

    GROUPS_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    candidates: list[dict] = []

    for accession in wanted:
        meta = registry.loc[accession]
        if accession in NO_CONTRAST:
            candidates.append({
                "accession": accession, "species": meta["species"],
                "n_samples": meta["n_samples"], "year": meta["year"],
                "assay": meta["assay"], "verdict": "exclude",
                "reason_code": NO_CONTRAST[accession],
                "reason_source": "single_cell_contrast_audit",
                "n_case": "", "n_control": "", "title": meta["title"]})
            logger.info("%-10s exclude  %s", accession, NO_CONTRAST[accession])
            continue

        rule = CONTRASTS.get(accession)
        if rule is None:
            logger.warning("%-10s no contrast rule and not in NO_CONTRAST", accession)
            continue

        rows, status = assign(accession, rule)
        if status != "assigned":
            candidates.append({
                "accession": accession, "species": meta["species"],
                "n_samples": meta["n_samples"], "year": meta["year"],
                "assay": meta["assay"], "verdict": "exclude",
                "reason_code": status, "reason_source": "step18_group_assignment",
                "n_case": "", "n_control": "", "title": meta["title"]})
            logger.error("%-10s %s", accession, status)
            continue

        all_rows.extend(rows)
        verdict, reason, n_case, n_control = study_verdict(rows)
        candidates.append({
            "accession": accession, "species": meta["species"],
            "n_samples": meta["n_samples"], "year": meta["year"],
            "assay": meta["assay"], "verdict": verdict, "reason_code": reason,
            "reason_source": "step18_group_assignment",
            "n_case": n_case, "n_control": n_control, "title": meta["title"]})

        placed = [r for r in rows if r["group"]]
        subjects = {r["subject"] for r in placed}
        merged = len(placed) - len(subjects)
        logger.info("%-10s %-7s %2d case / %2d control from %d GSM%s",
                    accession, verdict, n_case, n_control, len(placed),
                    f" ({merged} summed into a repeated subject)" if merged else "")

        if verdict == "include":
            frame = pd.DataFrame(placed)[["sample_id", "group", "subject", "title"]]
            frame.to_csv(GROUPS_DIR / f"{accession}_groups.csv", index=False)

    if all_rows:
        record = PRISMA / "single_cell_group_assignments.csv"
        pd.DataFrame(all_rows).to_csv(record, index=False)
        logger.info("%d GSM decisions -> %s", len(all_rows), record)

    if candidates:
        sheet = PRISMA / "single_cell_candidates.csv"
        pd.DataFrame(candidates).to_csv(sheet, index=False)
        included = [c for c in candidates if c["verdict"] == "include"]
        logger.info("%d of %d studies included -> %s",
                    len(included), len(candidates), sheet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
