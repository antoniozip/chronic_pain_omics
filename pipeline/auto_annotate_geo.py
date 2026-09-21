"""Automatically generate case/control sample annotations from GEO series matrix files.

Reads `*_series_matrix.txt.gz` files in `data/raw/geo_cache/` and emits
`{accession}_groups.csv` (sample_id, group) files where group ∈ {"case","control",""}.

Strategy:
  1. Parse !Sample_geo_accession, !Sample_title, !Sample_source_name_ch1,
     and all !Sample_characteristics_ch1 rows from the matrix header.
  2. For each candidate "signal column" (characteristics with >1 unique value,
     plus the title field), score each sample value against:
        - PAIN_TOKENS  → case
        - CTRL_TOKENS  → control
  3. Pick the column that produces the cleanest case/control split:
        - both classes have ≥ MIN_PER_GROUP samples,
        - the split covers the highest fraction of total samples,
        - tie-break: characteristics row preferred over title.
  4. Mark non-matching samples as "" (excluded).
  5. Skip when no column produces a valid split (genotype-only, single class,
     dose-response without ctrl, etc.) — caller will need manual annotation.

Outputs:
  data/raw/geo_cache/{accession}_groups.csv          for each annotated study
  data/raw/geo_cache/auto_annotation_log.csv         summary log
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------

# Tokens treated as CASE (pain / disease / injury). Order matters for logging;
# matching is whole-word or punctuation-bounded where helpful.
PAIN_TOKENS = (
    # Neuropathic models
    "sni", "cci", "snl", "psnl", "spnl", "spnt", "ion-cci", "ioncci",
    "spared nerve injury", "chronic constriction", "nerve ligation",
    "nerve injury", "neuropathic", "neuropathy", "axotomy",
    # Inflammatory models
    "cfa", "carrageenan", "carragenan", "formalin", "capsaicin",
    "complete freund", "freund's adjuvant", "monoiodoacetate", "mia",
    # Disease / clinical phenotypes
    "fibromyalgia", "osteoarthritis", "rheumatoid", "ibs", "endometriosis",
    "diabetic", "paclitaxel", "oxaliplatin", "vincristine", "chemotherapy",
    "cipn", "bone cancer", "cancer pain", "tumour", "tumor pain",
    "post-herpetic", "postherpetic", "phn", "hiv-sn",
    "crps", "complex regional",
    # Condition names the 2026-08-28 search amendment added to the query and
    # never to this table, so the retrieval could find them and the annotation
    # could not.
    "low back pain", "sciatica", "temporomandibular",
    "interstitial cystitis", "painful bladder", "bladder pain",
    # "hunner" is deliberately absent. It names a *lesion type* inside
    # interstitial cystitis, not a diagnosis, and with the negation rule it
    # manufactured a 25 v 25 contrast in GSE238208 out of Hunner-lesion and
    # non-lesional biopsies taken from the same bladders.
    # Generic pain wording
    "chronic pain", "pain patient", "pain group", "pain model",
    "allodynia", "hyperalgesia", "nociception", "nociceptive",
    "neuralgia", "trigeminal",
    # Surgery / injury (treated as case when paired with sham/naive)
    "ipsilateral", "contralateral",
    "surgery", "ligated", "injured", "lesion",
)

# Tokens treated as CONTROL
CTRL_TOKENS = (
    "sham", "naive", "naïve", "naieve",
    "control", "ctrl", "ctr", "ctl",
    " cont ", "_cont", "cont_", "(cont)",
    "healthy", "untreated", "vehicle", "saline", "pbs",
    "baseline", "normal", "no treatment", "no-treatment", "no_treatment",
    "uninjured", "no injury", "no surgery",
    "mock", "wild type", "wild-type", "wildtype",
)

# Tokens that, if present in the column name (key), make us skip the column —
# the variation is technical/demographic, not disease.
SKIP_KEY_TOKENS = (
    "tissue", "cell type", "cell line", "celltype",
    # Anatomical, not diagnostic. GSE238208 biopsies Hunner lesions and
    # non-lesional mucosa from the same bladders: "biopsied site" separates
    # perfectly, and reading it as a contrast makes a within-patient site
    # comparison look like disease against control.
    "site", "sites", "biopsy", "biopsied", "location", "region",
    "anatomic", "anatomical",
    "sex", "gender", "age", "strain", "background",
    "platform", "library", "kit", "instrument",
    "batch", "replicate", "rep ", "biological replicate",
    "cohort", "donor", "subject id", "sample id", "patient id",
    "fraction", "spike", "dropseq day", "scrna", "passage",
    "rin", "mapping",
)

# Minimum per-group sample count required to keep a study auto-annotated.
# Set to 2: downstream pipeline (cfg.da.min_samples_per_group) does final filtering.
MIN_PER_GROUP = 2

# Tokens that, when present in the column KEY, indicate the column splits a
# pain phenotype from a baseline. Any non-control value in such a column is
# treated as a case sample (key-driven inference).
PAIN_KEY_TOKENS = (
    "pain", "injury", "lesion", "surgery", "neuropath",
    "post sni", "post cci", "post snl", "days post",
    "model", "condition", "disease",
)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_row(line: str) -> tuple[str, list[str]]:
    """Split a !Sample_xxx tab-separated row, stripping quotes around values."""
    parts = line.rstrip("\n").split("\t")
    key = parts[0].lstrip("!")
    values = [p.strip().strip('"') for p in parts[1:]]
    return key, values


@dataclass
class SeriesMatrix:
    accession: str          # GSEnnnnn (without -GPLx)
    file_name: str
    gsm_ids: list[str]
    titles: list[str]
    source_names: list[str]
    characteristics: list[tuple[str, list[str]]]  # [(key, [values per sample])]

    @property
    def n(self) -> int:
        return len(self.gsm_ids)


def parse_series_matrix(path: Path) -> SeriesMatrix | None:
    """Parse the header of a GEO series_matrix.txt.gz file."""
    gsm_ids: list[str] = []
    titles: list[str] = []
    source_names: list[str] = []
    raw_chars: list[list[str]] = []  # list of [values per sample]

    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("!series_matrix_table_begin"):
                    break
                if not line.startswith("!Sample_"):
                    continue
                key, values = _parse_row(line)
                if key == "Sample_geo_accession":
                    gsm_ids = values
                elif key == "Sample_title":
                    titles = values
                elif key == "Sample_source_name_ch1":
                    source_names = values
                elif key.startswith("Sample_characteristics_ch1"):
                    raw_chars.append(values)
    except OSError as exc:
        logger.warning("Could not read %s: %s", path.name, exc)
        return None

    if not gsm_ids:
        return None

    # Convert each characteristics row "key: value" → ("key", ["v1","v2",...]).
    chars: list[tuple[str, list[str]]] = []
    for row in raw_chars:
        if not row:
            continue
        key_guess = _infer_key(row)
        # Strip the "key: " prefix from values to keep just the variable part
        clean = [_strip_key_prefix(v, key_guess) for v in row]
        chars.append((key_guess, clean))

    # Accession from filename, e.g. GSE111216_series_matrix.txt.gz → GSE111216
    stem = path.name.replace(".txt.gz", "")
    accession = re.sub(r"-GPL\w+$", "", stem.replace("_series_matrix", ""))

    if not source_names:
        source_names = [""] * len(gsm_ids)
    if not titles:
        titles = [""] * len(gsm_ids)

    return SeriesMatrix(
        accession=accession,
        file_name=path.name,
        gsm_ids=gsm_ids,
        titles=titles,
        source_names=source_names,
        characteristics=chars,
    )


def _infer_key(values: list[str]) -> str:
    """Most GEO characteristics values look like 'key: value'; return the key."""
    keys = []
    for v in values:
        if ":" in v:
            keys.append(v.split(":", 1)[0].strip().lower())
    if not keys:
        return "characteristic"
    # Pick the most common key
    return max(set(keys), key=keys.count)


def _strip_key_prefix(value: str, key: str) -> str:
    """Drop a leading 'key:' from a value if present."""
    if ":" in value:
        head, _, tail = value.partition(":")
        if head.strip().lower() == key.strip().lower():
            return tail.strip()
    return value.strip()


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

# Words that negate a diagnosis. A control arm routinely names the disease it
# does not have -- "no endometriosis", "non-endometriosis", "disease-free" --
# and the disease token then matched, so *both* arms classified as case and the
# study reported no split at all. That single defect hid a body of
# endometriosis case/control studies behind `no_case_control_split_in_series_matrix`.
_NEG_PREFIX = r"(?:non|no|without|w/o|free\s+of|negative\s+for|absence\s+of)"
#: Generic disease words, for negations that name no specific condition.
_GENERIC_DISEASE = ("disease", "diseased", "pain", "lesion", "lesions",
                    "symptom", "symptoms", "pathology")


def _negates_disease(padded_lower: str) -> bool:
    """True when the value negates a disease term rather than carrying one."""
    for tok in tuple(PAIN_TOKENS) + _GENERIC_DISEASE:
        t = re.escape(tok.lower())
        if re.search(rf"(?<![a-z0-9]){_NEG_PREFIX}[-_ ]?{t}(?![a-z0-9])",
                     padded_lower):
            return True
        if re.search(rf"(?<![a-z0-9]){t}[-_ ]?free(?![a-z0-9])", padded_lower):
            return True
    return False


def _classify_value(val: str) -> str | None:
    """Return 'case', 'control', or None for a single sample value."""
    lo = " " + val.lower() + " "
    # Checked first: a negated disease term also matches that disease term.
    if _negates_disease(lo):
        return "control"
    has_pain = any(_token_present(tok, lo) for tok in PAIN_TOKENS)
    has_ctrl = any(_token_present(tok, lo) for tok in CTRL_TOKENS)
    if has_pain and not has_ctrl:
        return "case"
    if has_ctrl and not has_pain:
        return "control"
    return None


def _token_present(tok: str, padded_lower: str) -> bool:
    """Substring match; phrases match anywhere, short codes need word boundary."""
    tok = tok.lower()
    if " " in tok or "-" in tok or "'" in tok:
        return tok in padded_lower
    # Short alphanumeric token: require non-letter neighbours so 'sni' doesn't
    # match 'thoracoscopic' etc.
    pattern = r"(?<![a-z0-9])" + re.escape(tok) + r"(?![a-z0-9])"
    return re.search(pattern, padded_lower) is not None


def _column_should_skip(key: str) -> bool:
    """True when the characteristics key is purely technical/demographic.

    Matching follows `_token_present`: phrases anywhere, single tokens on a
    boundary. A plain substring test skipped "disease stage" because "age" is
    inside "stage", which discarded the disease column of GSE141549 -- 408
    samples, the largest study in the arm.
    """
    k = " " + key.lower() + " "
    return any(_token_present(t, k) for t in SKIP_KEY_TOKENS)


@dataclass
class Candidate:
    source: str                      # "char:<key>" or "title" or "source_name"
    labels: list[str]                # one of "case"/"control"/"" per sample
    n_case: int
    n_ctrl: int
    coverage: float                  # (n_case+n_ctrl)/n


def _score(labels: list[str]) -> tuple[int, int, float]:
    n_case = sum(1 for x in labels if x == "case")
    n_ctrl = sum(1 for x in labels if x == "control")
    n = len(labels)
    return n_case, n_ctrl, (n_case + n_ctrl) / n if n else 0.0


def _key_is_pain_signal(key: str) -> bool:
    k = " " + key.lower() + " "
    return any(t in k for t in PAIN_KEY_TOKENS)


def _candidate_from_values(source: str, values: list[str],
                           key_drives_case: bool = False) -> Candidate:
    """Classify each value. When key_drives_case=True, non-control non-empty
    values default to 'case' instead of None."""
    labels: list[str] = []
    for v in values:
        cls = _classify_value(v)
        if cls is None and key_drives_case and v.strip():
            # The column key itself signals the pain axis; anything that isn't
            # explicitly control is presumed to be case.
            labels.append("case")
        else:
            labels.append(cls or "")
    n_case, n_ctrl, cov = _score(labels)
    return Candidate(source=source, labels=labels, n_case=n_case,
                     n_ctrl=n_ctrl, coverage=cov)


def annotate(sm: SeriesMatrix) -> Candidate | None:
    """Return the best (source, label vector) for a series matrix, or None."""
    candidates: list[Candidate] = []

    # Characteristics rows
    for key, vals in sm.characteristics:
        if _column_should_skip(key):
            continue
        # Need real variation across samples
        if len(set(v.lower() for v in vals)) < 2:
            continue
        key_drives_case = _key_is_pain_signal(key)
        candidates.append(_candidate_from_values(f"char:{key}", vals,
                                                 key_drives_case=key_drives_case))

    # Source name & title as fallbacks
    if len(set(sm.source_names)) >= 2:
        candidates.append(_candidate_from_values("source_name", sm.source_names))
    if len(set(sm.titles)) >= 2:
        candidates.append(_candidate_from_values("title", sm.titles))

    # Keep only candidates with both classes meeting MIN_PER_GROUP
    viable = [c for c in candidates
              if c.n_case >= MIN_PER_GROUP and c.n_ctrl >= MIN_PER_GROUP]
    if not viable:
        return None

    # Prefer characteristics rows over title; then highest coverage; then
    # smaller class-size imbalance.
    def rank(c: Candidate) -> tuple:
        source_pref = 0 if c.source.startswith("char:") else (1 if c.source == "source_name" else 2)
        balance = -abs(c.n_case - c.n_ctrl)
        return (source_pref, -c.coverage, -balance)

    viable.sort(key=rank)
    return viable[0]


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def write_groups_csv(out_path: Path, sm: SeriesMatrix, cand: Candidate) -> None:
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "group"])
        for gsm, label in zip(sm.gsm_ids, cand.labels):
            w.writerow([gsm, label])


SINGLE_CELL_STUDIES = (Path(__file__).resolve().parent.parent
                       / "conf" / "analysis" / "single_cell_studies.csv")
DRY_RUN_LOG = (Path(__file__).resolve().parent.parent
               / "temp" / "auto_annotation_dry_run.csv")


def single_cell_accessions() -> set[str]:
    """Accessions belonging to the single-cell arm.

    Every one of them is *also* in the transcriptomic screening record, and
    `{accession}_groups.csv` in this shared cache is exactly what
    `run_transcriptomics_da()` reads. Writing one here wires the two arms
    together the moment anyone flips a verdict, which is why the arm keeps its
    own groups under data/interim/single_cell/groups/.
    """
    if not SINGLE_CELL_STUDIES.exists():
        return set()
    with open(SINGLE_CELL_STUDIES, newline="") as fh:
        return {row["accession"] for row in csv.DictReader(fh)
                if row.get("accession")}


def run(geo_cache_dir: Path, overwrite: bool = False,
        dry_run: bool = False) -> None:
    matrices = sorted(geo_cache_dir.glob("*_series_matrix.txt.gz"))
    reserved = single_cell_accessions()
    log_rows: list[dict] = []
    n_written = 0
    n_skipped = 0
    n_existing = 0

    # Group multi-platform variants by accession; pick the file with the most samples
    by_accession: dict[str, list[Path]] = {}
    for path in matrices:
        stem = path.name.replace("_series_matrix.txt.gz", "")
        acc = re.sub(r"-GPL\w+$", "", stem)
        by_accession.setdefault(acc, []).append(path)

    for accession, paths in by_accession.items():
        out_path = geo_cache_dir / f"{accession}_groups.csv"
        if accession in reserved:
            n_skipped += 1
            log_rows.append({
                "accession": accession, "status": "skipped_single_cell_arm",
                "source": "", "n_samples": "", "n_case": "", "n_control": "",
                "reason": "belongs to the single-cell arm; a groups file here "
                          "would wire the two arms together",
            })
            continue
        if out_path.exists() and not overwrite:
            n_existing += 1
            log_rows.append({
                "accession": accession, "status": "skipped_existing",
                "source": "", "n_samples": "", "n_case": "", "n_control": "",
                "reason": "manual or prior auto annotation present",
            })
            continue

        # Pick the largest matrix for this accession
        parsed = [p for p in (parse_series_matrix(pp) for pp in paths) if p]
        if not parsed:
            n_skipped += 1
            log_rows.append({
                "accession": accession, "status": "skipped_parse_error",
                "source": "", "n_samples": "", "n_case": "", "n_control": "",
                "reason": "could not parse series matrix",
            })
            continue
        sm = max(parsed, key=lambda s: s.n)

        cand = annotate(sm)
        if cand is None:
            n_skipped += 1
            log_rows.append({
                "accession": accession, "status": "no_split",
                "source": "", "n_samples": sm.n, "n_case": "", "n_control": "",
                "reason": "no candidate column produced a case/control split",
            })
            continue

        if not dry_run:
            write_groups_csv(out_path, sm, cand)
        n_written += 1
        log_rows.append({
            "accession": accession,
            "status": "would_annotate" if dry_run else "annotated",
            "source": cand.source, "n_samples": sm.n,
            "n_case": cand.n_case, "n_control": cand.n_ctrl,
            "reason": "",
        })

    # Write the summary log. A dry run keeps it out of the cache too, so an
    # audit leaves the tree exactly as it found it. The dry-run path is fixed
    # to the repository rather than derived from the cache directory, so it
    # does not land somewhere arbitrary when the cache is a temporary one.
    log_path = DRY_RUN_LOG if dry_run else geo_cache_dir / "auto_annotation_log.csv"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "accession", "status", "source", "n_samples",
            "n_case", "n_control", "reason",
        ])
        w.writeheader()
        w.writerows(log_rows)

    print("Auto-annotation complete:" + ("  (DRY RUN, nothing written)"
                                          if dry_run else ""))
    print(f"  {'Would write' if dry_run else 'Written'}:  {n_written}")
    print(f"  Skipped:  {n_skipped}")
    print(f"  Existing: {n_existing}")
    print(f"  Log:      {log_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geo-cache-dir", type=Path,
                        default=Path("data/raw/geo_cache"))
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing *_groups.csv files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be annotated, write nothing")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run(args.geo_cache_dir, overwrite=args.overwrite, dry_run=args.dry_run)
