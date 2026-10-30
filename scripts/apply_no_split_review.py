"""Apply the hand review of the no-split exclusions to the screening sheet.

`review_no_split.py` records one decision per series in
`literature/prisma/no_split_review.csv` and never writes the sheet. This is the
deliberate second step. For each `contrast_exists` series it

  * recomputes the case/control labels from the saved rule, refusing a series
    whose counts no longer match what the reviewer saw,
  * writes `data/raw/geo_cache/{acc}_groups.csv` (case and control samples only,
    so every other sample is an exclusion rather than a silent drop), and
  * moves the sheet row from `exclude` to `pending / usable_pending_file_check`,
    the state "cleared on design, files unchecked". It is not `include`: only a
    step 05 run that writes an effect table earns that.

`no_contrast` rows are left as they were. Dry-run unless given `--apply`; an
existing groups file is never overwritten.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import review_no_split as rv  # noqa: E402

logger = logging.getLogger(__name__)

SOURCE = "no_split_review_2026_10_21"
PENDING = "usable_pending_file_check"


def groups_for(root: Path, row: dict[str, str]) -> list[tuple[str, str]]:
    """Case/control GSMs for one saved `contrast_exists` review, counts re-checked."""
    rule = json.loads(row["rule_json"])
    detail = rv._load(str(root), row["accession"])
    fields = {f["id"]: f for f in detail["fields"]}
    labels = rv.apply_rule(fields[rule["field_id"]]["values"], rule)
    got = (labels["n_case"], labels["n_control"])
    saved = (int(row["n_case"]), int(row["n_control"]))
    if got != saved or labels["n_conflict"]:
        raise ValueError(f"{row['accession']}: rule gives {got} (conflicts "
                         f"{labels['n_conflict']}), review saved {saved}")
    return [(g, lab) for g, lab in zip(detail["gsm"], labels["labels"], strict=True)
            if lab in ("case", "control")]


def apply(root: Path, reviews: Path, write: bool) -> list[str]:
    sheet = root / rv.SHEET
    raw = sheet.read_bytes().decode("utf-8")
    rows = list(csv.DictReader(raw.splitlines()))
    fields = list(rows[0])
    by_acc = {r["accession"]: r for r in rows}
    done: list[str] = []
    for acc, rev in sorted(rv.load_reviews(reviews).items()):
        if rev["decision"] != "contrast_exists":
            continue
        row = by_acc[acc]
        if row["reason_code"] != rv.REASON:
            raise ValueError(f"{acc}: sheet reason is {row['reason_code']!r}, not the "
                             "no-split exclusion; refusing to overwrite another decision")
        out = root / rv.CACHE / f"{acc}_groups.csv"
        if out.exists():
            raise FileExistsError(f"{out} exists; remove it deliberately to re-derive")
        pairs = groups_for(root, rev)
        row.update(verdict="pending", reason_code=PENDING, reason_source=SOURCE,
                   annotation_status="annotated", n_case=rev["n_case"],
                   n_control=rev["n_control"])
        if write:
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["sample_id", "group"])
                w.writerows(pairs)
        done.append(acc)
    if write and done:
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\r\n")
        w.writeheader()
        w.writerows(rows)
        sheet.write_bytes(buf.getvalue().encode("utf-8"))
    return done


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=REPO_ROOT)
    p.add_argument("--reviews", type=Path, default=None)
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()
    done = apply(a.root, a.reviews or a.root / rv.REVIEWS, a.apply)
    logger.info("%s %d series to pending/%s: %s", "moved" if a.apply else "would move",
                len(done), PENDING, " ".join(done))


if __name__ == "__main__":
    main()
