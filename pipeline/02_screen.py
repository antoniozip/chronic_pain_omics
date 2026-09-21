"""Screening pipeline step.

Reads all JSONL hit files produced by 01_search_literature.py, deduplicates
across runs, and produces a screening sheet for title/abstract review.

Screening decisions are written back into the same CSV so the file is the
single source of truth for the PRISMA flow.

Usage:
    uv run python pipeline/02_screen.py                    # generate screening sheet
    uv run python pipeline/02_screen.py --status           # print PRISMA counts
    uv run python pipeline/02_screen.py --apply decisions.csv  # merge coded sheet back
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Literal

from cp_multiomics.vocabulary import PainVocabulary, load_pain_vocabulary

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
RAW_LITERATURE_DIR = REPO_ROOT / "data" / "raw" / "literature"
SCREENING_CSV = REPO_ROOT / "literature" / "prisma" / "screening.csv"

Decision = Literal["include", "exclude", "uncertain", ""]

EXCLUSION_REASONS = [
    "no_omics_data",
    "wrong_condition",          # not chronic pain
    "wrong_species",            # not human or target animal model
    "review_no_primary_data",
    "conference_abstract_only",
    "retracted",
    "duplicate",
    "no_omics_term",            # auto: no omics vocabulary term in title/abstract
    "not_pain_phenotype",       # auto: no pain vocabulary term in title/abstract
    "other",
]

# Reason codes only ever written by apply_auto_exclusions(). A row is safe to
# reset only if it still bears BOTH the "auto" tag AND one of these reasons —
# a human who overrode an auto row (e.g. to "include") via --apply may leave
# the stale screened_by="auto" tag in place, and must not be silently wiped.
AUTO_REASONS = {"no_omics_term", "not_pain_phenotype"}

# Omics terms used by auto-screening. Deliberately broad: a false negative here
# silently discards a study, whereas a false positive only costs a human review.
# Covers: transcriptomics, genomics/GWAS, proteomics, metabolomics/lipidomics,
# epigenomics (ATAC-seq/ChIP-seq/CUT&RUN/CUT&Tag, methylation, bisulfite),
# exome/genome sequencing, ribosome profiling, and single-cell assays.
# Bare 3-char acronyms "wes"/"wgs" are deliberately NOT included as raw
# substrings: "wes" collides with common words (e.g. "western"), which would
# flood human review with noise. Coverage instead comes from "exome"
# (>=5 chars, safe as a substring) and the spelled-out "whole-exome"/
# "whole-genome" forms.
DEFAULT_OMICS_TERMS = (
    "atac-seq", "atac seq", "scatac",
    "bisulfite",
    "chip-seq", "chip seq",
    "cut&run", "cut&tag", "cut and run", "cut and tag",
    "epigenom", "epigenetic", "epimutation",
    "exome", "whole-exome", "whole exome", "whole-genome", "whole genome",
    "expression profiling",
    "genome-wide association", "gwas",
    "lipidom",
    "mass spectrometry", "metabolom", "methylat", "microarray",
    "proteom",
    "ribosome profiling",
    "rna-seq", "rnaseq",
    "sequencing", "single-cell", "single cell",
    "snp",
    "transcriptom",
)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class ScreeningRecord:
    uid: str
    source: str
    doi: str
    title: str
    authors: str
    journal: str
    year: str
    abstract: str
    query_label: str
    retrieved_at: str
    # screening fields (filled by reviewer)
    screening_decision: Decision = ""
    exclusion_reason: str = ""
    screened_by: str = ""
    screening_notes: str = ""

    @classmethod
    def from_hit(cls, hit: dict) -> ScreeningRecord:
        return cls(
            uid=hit.get("uid", ""),
            source=hit.get("source", ""),
            doi=hit.get("doi", ""),
            title=hit.get("title", ""),
            authors=hit.get("authors", ""),
            journal=hit.get("journal", ""),
            year=hit.get("year", ""),
            abstract=hit.get("abstract", ""),
            query_label=hit.get("query_label", ""),
            retrieved_at=hit.get("retrieved_at", ""),
        )

    def dedup_key(self) -> str:
        if self.doi:
            return f"doi:{self.doi.lower().strip()}"
        return f"{self.source}:{self.uid}"


# ---------------------------------------------------------------------------
# Load raw hits from all JSONL files
# ---------------------------------------------------------------------------

def load_all_hits(raw_dir: Path) -> list[ScreeningRecord]:
    records: list[ScreeningRecord] = []
    jsonl_files = sorted(raw_dir.glob("*.jsonl"))
    if not jsonl_files:
        logger.warning("No JSONL hit files found in %s", raw_dir)
        return records

    for path in jsonl_files:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    records.append(ScreeningRecord.from_hit(json.loads(line)))
        logger.info("Loaded %s (%d records so far)", path.name, len(records))

    return records


def deduplicate(records: list[ScreeningRecord]) -> tuple[list[ScreeningRecord], int]:
    seen: set[str] = set()
    unique: list[ScreeningRecord] = []
    for r in records:
        key = r.dedup_key()
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique, len(records) - len(unique)


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

_FIELDNAMES = [f.name for f in fields(ScreeningRecord)]


def load_existing_screening(csv_path: Path) -> dict[str, ScreeningRecord]:
    """Return existing screening records keyed by dedup_key."""
    if not csv_path.exists():
        return {}
    existing: dict[str, ScreeningRecord] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            r = ScreeningRecord(**{k: row.get(k, "") for k in _FIELDNAMES})
            existing[r.dedup_key()] = r
    return existing


def write_screening_csv(records: list[ScreeningRecord], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)
        writer.writeheader()
        for r in records:
            writer.writerow({f.name: getattr(r, f.name) for f in fields(r)})
    logger.info("Screening sheet written → %s (%d records)", csv_path, len(records))


# ---------------------------------------------------------------------------
# Merge coded decisions back
# ---------------------------------------------------------------------------

def apply_decisions(screening_csv: Path, decisions_csv: Path) -> None:
    """Merge screener-filled decisions CSV back into the master screening sheet.

    The decisions CSV must contain at minimum: uid, source, doi,
    screening_decision, exclusion_reason, screened_by, screening_notes.
    Unrecognised exclusion reasons are flagged as warnings.
    """
    existing = load_existing_screening(screening_csv)

    updated = 0
    with open(decisions_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            tmp = ScreeningRecord.from_hit(row)
            key = tmp.dedup_key()
            if key not in existing:
                logger.warning("Unknown record in decisions file (key=%s), skipping", key)
                continue
            decision = row.get("screening_decision", "").strip().lower()
            reason = row.get("exclusion_reason", "").strip()
            if decision not in ("include", "exclude", "uncertain", ""):
                logger.warning("Invalid decision '%s' for %s — set to uncertain", decision, key)
                decision = "uncertain"
            if reason and reason not in EXCLUSION_REASONS:
                logger.warning("Unknown exclusion reason '%s' for %s", reason, key)
            rec = existing[key]
            rec.screening_decision = decision  # type: ignore[assignment]
            rec.exclusion_reason = reason
            screened_by = row.get("screened_by", "")
            if screened_by == "auto":
                logger.warning(
                    "Reviewer-supplied row %s carries the reserved screened_by='auto' "
                    "sentinel — likely a stale round-tripped value on a human edit "
                    "(decision=%r). Value kept as-is; verify with the reviewer.",
                    key, decision,
                )
            rec.screened_by = screened_by
            rec.screening_notes = row.get("screening_notes", "")
            updated += 1

    write_screening_csv(list(existing.values()), screening_csv)
    logger.info("Applied %d decisions from %s", updated, decisions_csv)


# ---------------------------------------------------------------------------
# PRISMA status report
# ---------------------------------------------------------------------------

def print_status(screening_csv: Path) -> None:
    if not screening_csv.exists():
        logger.info("No screening sheet found at %s", screening_csv)
        return

    records = list(load_existing_screening(screening_csv).values())
    total = len(records)
    included = sum(1 for r in records if r.screening_decision == "include")
    excluded = sum(1 for r in records if r.screening_decision == "exclude")
    uncertain = sum(1 for r in records if r.screening_decision == "uncertain")
    pending = sum(1 for r in records if r.screening_decision == "")

    reason_counts: dict[str, int] = {}
    for r in records:
        if r.exclusion_reason:
            reason_counts[r.exclusion_reason] = reason_counts.get(r.exclusion_reason, 0) + 1

    print("\n=== PRISMA Screening Status ===")
    print(f"  Total records  : {total:>6}")
    print(f"  Included       : {included:>6}")
    print(f"  Excluded       : {excluded:>6}")
    print(f"  Uncertain      : {uncertain:>6}")
    print(f"  Pending review : {pending:>6}")
    if reason_counts:
        print("\n  Exclusion reasons:")
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"    {reason:<30} {count:>5}")
    print()


# ---------------------------------------------------------------------------
# Auto-screening (exclusion only — inclusion is always a human decision)
# ---------------------------------------------------------------------------

def auto_screen(
    record: ScreeningRecord,
    vocab: PainVocabulary,
    omics_terms: tuple[str, ...] = DEFAULT_OMICS_TERMS,
) -> tuple[str, str]:
    """Decide whether a record can be excluded without human review.

    Returns ("", "") when the record must be reviewed by a human. Never returns
    "include": PRISMA inclusion is always a human decision.

    Args:
        record: The screening record.
        vocab: The pain vocabulary.
        omics_terms: Lowercase omics substrings.

    Returns:
        (decision, exclusion_reason).
    """
    blob = f"{record.title} {record.abstract}".lower()

    if not vocab.matches(blob):
        return "exclude", "not_pain_phenotype"
    if not any(term in blob for term in omics_terms):
        return "exclude", "no_omics_term"
    return "", ""


def apply_auto_exclusions(
    screening_csv: Path,
    vocab: PainVocabulary,
    omics_terms: tuple[str, ...] = DEFAULT_OMICS_TERMS,
) -> dict[str, int]:
    """Fill blank decisions with auto-exclusions. Never overwrites a human decision."""
    existing = load_existing_screening(screening_csv)
    counts = {"excluded": 0, "left_for_human": 0, "skipped_already_decided": 0}

    for rec in existing.values():
        if rec.screening_decision:
            counts["skipped_already_decided"] += 1
            continue
        decision, reason = auto_screen(rec, vocab, omics_terms)
        if decision == "exclude":
            rec.screening_decision = decision  # type: ignore[assignment]
            rec.exclusion_reason = reason
            rec.screened_by = "auto"
            counts["excluded"] += 1
        else:
            counts["left_for_human"] += 1

    write_screening_csv(list(existing.values()), screening_csv)
    logger.info(
        "Auto-exclusion: %d excluded, %d left for human review, %d already decided",
        counts["excluded"], counts["left_for_human"], counts["skipped_already_decided"],
    )
    return counts


def reset_auto_exclusions(screening_csv: Path) -> int:
    """Clear every auto-made decision so it can be re-screened with an updated vocabulary.

    Resets screening_decision, exclusion_reason, and screened_by back to "" for
    every row that bears BOTH the "auto" tag AND an auto-only reason code
    (AUTO_REASONS) with decision "exclude" — i.e. a row that can only have
    been written by apply_auto_exclusions(). Gating on screened_by alone is
    not safe: --apply round-trips screened_by verbatim, so a human who
    overrides an auto-excluded row (e.g. to "include") without also editing
    screened_by leaves the stale "auto" tag in place. Requiring the reason
    code too means a human edit to screening_decision or exclusion_reason
    takes the row out of scope for reset, even if the tag was never touched.

    Args:
        screening_csv: Path to the master screening sheet.

    Returns:
        The number of rows reset.
    """
    existing = load_existing_screening(screening_csv)
    reset_count = 0
    for rec in existing.values():
        if rec.screened_by != "auto":
            continue
        if rec.screening_decision != "exclude" or rec.exclusion_reason not in AUTO_REASONS:
            logger.info(
                "Preserving row with screened_by='auto' but decision=%r reason=%r "
                "(looks like a human override that left the stale 'auto' tag) — "
                "not reset",
                rec.screening_decision, rec.exclusion_reason,
            )
            continue
        rec.screening_decision = ""  # type: ignore[assignment]
        rec.exclusion_reason = ""
        rec.screened_by = ""
        reset_count += 1

    write_screening_csv(list(existing.values()), screening_csv)
    logger.info("Reset %d auto-made decisions back to pending", reset_count)
    return reset_count


def audit_sample(screening_csv: Path, n: int = 200, seed: int = 42) -> list[ScreeningRecord]:
    """Return a deterministic random sample of auto-excluded records for human audit.

    The false-negative rate estimated from this sample MUST be reported in the
    manuscript alongside the auto-exclusion count.
    """
    import random

    records = [
        r for r in load_existing_screening(screening_csv).values() if r.screened_by == "auto"
    ]
    records.sort(key=lambda r: r.uid)
    rng = random.Random(seed)
    return rng.sample(records, min(n, len(records)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--raw-dir",
        type=Path,
        default=RAW_LITERATURE_DIR,
        help="Directory containing JSONL hit files from step 01",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=SCREENING_CSV,
        help="Path to write (or update) the screening CSV",
    )
    p.add_argument(
        "--apply",
        type=Path,
        default=None,
        metavar="DECISIONS_CSV",
        help="Merge a reviewer-filled decisions CSV back into the master sheet",
    )
    p.add_argument(
        "--status",
        action="store_true",
        help="Print PRISMA screening counts and exit",
    )
    p.add_argument(
        "--auto-exclude",
        action="store_true",
        help="Fill blank decisions with vocabulary-based auto-exclusions (never auto-includes)",
    )
    p.add_argument(
        "--reset-auto",
        action="store_true",
        help="Clear all auto-made decisions (screened_by == 'auto') back to pending; "
             "never touches a human decision",
    )
    p.add_argument(
        "--audit-sample",
        type=int,
        default=0,
        metavar="N",
        help="Print N randomly sampled auto-excluded records for human audit",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.status:
        print_status(args.output)
        return

    if args.apply:
        apply_decisions(args.output, args.apply)
        print_status(args.output)
        return

    if args.auto_exclude:
        apply_auto_exclusions(args.output, load_pain_vocabulary())
        print_status(args.output)
        return

    if args.reset_auto:
        n = reset_auto_exclusions(args.output)
        print(f"Reset {n} auto-made decisions to pending.")
        print_status(args.output)
        return

    if args.audit_sample:
        for rec in audit_sample(args.output, n=args.audit_sample):
            print(f"{rec.uid}\t{rec.exclusion_reason}\t{rec.title[:90]}")
        return

    # Default: generate / refresh screening sheet from raw hits
    all_hits = load_all_hits(args.raw_dir)
    if not all_hits:
        logger.info("No hits to screen. Run 01_search_literature.py first.")
        return

    unique, n_dupes = deduplicate(all_hits)
    logger.info("Loaded %d hits, removed %d cross-run duplicates → %d unique",
                len(all_hits), n_dupes, len(unique))

    # Preserve existing decisions when refreshing
    existing = load_existing_screening(args.output)
    preserved = 0
    for r in unique:
        key = r.dedup_key()
        if key in existing:
            prev = existing[key]
            r.screening_decision = prev.screening_decision
            r.exclusion_reason = prev.exclusion_reason
            r.screened_by = prev.screened_by
            r.screening_notes = prev.screening_notes
            if prev.screening_decision:
                preserved += 1

    if preserved:
        logger.info("Preserved %d existing decisions", preserved)

    write_screening_csv(unique, args.output)
    print_status(args.output)


if __name__ == "__main__":
    main()
