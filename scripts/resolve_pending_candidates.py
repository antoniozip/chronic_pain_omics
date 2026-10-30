#!/usr/bin/env python
"""Apply screener decisions to pending transcriptomic candidate rows.

Accession-keyed sibling of `02_screen.py --apply`, for the PRISMA candidate
sheet the per-study layer reads: `literature/prisma/transcriptomics_candidates.csv`.

Only rows whose verdict is `pending` may change. Every exclude must carry a
reason code; an include must not. The sheet is written atomically (tmp +
os.replace), so a crash cannot leave it half-written.

Usage:
    python scripts/resolve_pending_candidates.py \
        --decisions literature/prisma/pending_resolution_A.csv \
        --sheet literature/prisma/transcriptomics_candidates.csv

The decisions CSV needs at least the columns: accession, verdict, reason_code.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SHEET = REPO_ROOT / "literature" / "prisma" / "transcriptomics_candidates.csv"

REQUIRED_DECISION_COLS = {"accession", "verdict"}
REQUIRED_SHEET_COLS = {"accession", "verdict", "reason_code"}


def _fail(message: str) -> None:
    print(message, file=sys.stderr)
    sys.exit(1)


def resolve(sheet: pd.DataFrame, decisions: pd.DataFrame,
            allow_include_to_exclude: bool = False) -> pd.DataFrame:
    """Return the sheet with the pending rows updated; raise ValueError on bad input."""
    missing = REQUIRED_SHEET_COLS - set(sheet.columns)
    if missing:
        raise ValueError(f"sheet lacks column(s): {sorted(missing)}")
    missing = REQUIRED_DECISION_COLS - set(decisions.columns)
    if missing:
        raise ValueError(f"decisions CSV lacks column(s): {sorted(missing)}")

    out = sheet.copy()
    changed = 0
    for _, row in decisions.iterrows():
        acc = str(row["accession"]).strip()
        verdict = str(row["verdict"]).strip().lower()
        raw_reason = row.get("reason_code")
        reason = "" if raw_reason is None or pd.isna(raw_reason) else str(raw_reason).strip()

        mask = out["accession"].astype(str) == acc
        if not mask.any():
            raise ValueError(f"unknown accession {acc!r}: not present in the sheet")
        current = out.loc[mask, "verdict"].iloc[0]
        if current == "include" and verdict == "exclude" and allow_include_to_exclude:
            pass  # provisional include from a pending-resolution pass; re-decided below
        elif current != "pending":
            raise ValueError(
                f"{acc}: current verdict is {current!r}, not pending; refusing to re-judge")
        if verdict not in ("include", "exclude"):
            raise ValueError(f"{acc}: invalid verdict {verdict!r}")
        if verdict == "exclude" and not reason:
            raise ValueError(f"{acc}: an exclude decision requires a reason_code")
        if verdict == "include" and reason:
            raise ValueError(f"{acc}: an include decision must not carry a reason_code")

        out.loc[mask, "verdict"] = verdict
        out.loc[mask, "reason_code"] = reason
        changed += 1

    if changed == 0:
        raise ValueError("no decision rows to apply")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--decisions", type=Path, required=True)
    ap.add_argument("--sheet", type=Path, default=DEFAULT_SHEET)
    ap.add_argument("--allow-include-to-exclude", action="store_true",
                    help="permit re-deciding an include row back to exclude "
                         "(provisional includes from a pending-resolution pass)")
    args = ap.parse_args()

    sheet = pd.read_csv(args.sheet)
    decisions = pd.read_csv(args.decisions)

    try:
        out = resolve(sheet, decisions,
                        allow_include_to_exclude=args.allow_include_to_exclude)
    except ValueError as exc:
        _fail(f"error: {exc}")

    before = sheet["verdict"].value_counts().to_dict()
    after = out["verdict"].value_counts().to_dict()
    tmp = args.sheet.with_suffix(args.sheet.suffix + ".tmp")
    out.to_csv(tmp, index=False)
    os.replace(tmp, args.sheet)

    print(f"applied {len(decisions)} decisions")
    print("verdicts before:", before)
    print("verdicts after: ", after)
    print(f"pending remaining: {after.get('pending', 0)}")


if __name__ == "__main__":
    main()
