"""Apply the hand mapping of count columns to the cache and the screening sheet.

`map_count_columns.py` records one decision per series in
`literature/prisma/count_column_mapping.csv` and writes nothing else. This is the
deliberate second step. For each `mapping_ok` series it

  * writes `data/raw/geo_cache/{acc}/{acc}_raw_count_merged.txt.gz`: the gene
    column and the paired sample columns, renamed to GSM ids, which is the one
    naming that `{acc}_groups.csv` can be matched against;
  * moves every other file of that folder, the deposit itself included, to
    `temp/count_mapping_originals/{acc}/`, because step 05 picks a count file by
    tier and has no preference by name, so a second candidate makes the choice
    alphabetical;
  * rewrites `{acc}_groups.csv` without the samples the reviewer dropped (the
    original is kept beside the others);
  * leaves the sheet row at `pending / usable_pending_file_check`, with the
    arm counts the reviewer's mapping gives. It is not `include`: only a step 05
    run that writes an effect table earns that.

An `ineligible` decision moves the sheet row to `exclude` with the reviewer's
reason code and touches no data. `unclear` changes nothing. Dry-run unless given
`--apply`; a series whose sheet row is no longer pending for this reason is
refused rather than overwritten.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import map_count_columns as mc  # noqa: E402

rv = mc.rv
logger = logging.getLogger(__name__)

SOURCE = "count_column_mapping_2026_10_23"
PENDING = "usable_pending_file_check"
ORIGINALS = Path("temp/count_mapping_originals")


def build_table(df: pd.DataFrame, rev: dict[str, str]) -> pd.DataFrame:
    """The reviewer's columns, renamed to GSM ids, on the scale step 05 reads.

    `log2` values are returned to the linear scale as 2**x - 1: the inverse of the
    log2(x + 1) the deposits use, which the zeros in them confirm. Duplicate gene
    labels keep the row with the highest mean, since summing is only right for
    counts and the choice should not depend on the scale.
    """
    mapping: dict[str, str] = json.loads(rev["mapping_json"])
    drop = set(json.loads(rev["drop_json"] or "[]"))
    kept = {c: g for c, g in mapping.items() if g not in drop}
    body = df[list(kept)].apply(pd.to_numeric, errors="coerce").rename(columns=kept)
    if rev["scale"] == "log2":
        body = (2.0 ** body - 1.0).clip(lower=0)
    genes = df[rev["gene_col"]].astype("string").str.strip()
    body.insert(0, "gene_id", genes)
    body = body[body["gene_id"].notna() & (body["gene_id"] != "")]
    order = body.drop(columns="gene_id").mean(axis=1).sort_values(ascending=False).index
    body = body.loc[order].drop_duplicates("gene_id").sort_index()
    return body.reset_index(drop=True)


def _write_sheet(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\r\n")
    w.writeheader()
    w.writerows(rows)
    path.write_bytes(buf.getvalue().encode("utf-8"))


def _move_aside(root: Path, accession: str) -> list[str]:
    folder = root / rv.CACHE / accession
    keep = root / ORIGINALS / accession
    moved = []
    for p in sorted(folder.glob("*")):
        if p.is_file():
            keep.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(keep / p.name))
            moved.append(p.name)
    return moved


def apply(root: Path, reviews: Path, write: bool) -> dict[str, list[str]]:
    sheet = root / rv.SHEET
    rows = list(csv.DictReader(sheet.read_bytes().decode("utf-8").splitlines()))
    fields = list(rows[0])
    by_acc = {r["accession"]: r for r in rows}
    done: dict[str, list[str]] = {"mapped": [], "excluded": []}
    for acc, rev in sorted(rv.load_reviews(reviews).items()):
        if rev["decision"] not in ("mapping_ok", "ineligible"):
            continue
        row = by_acc[acc]
        if row["reason_code"] != mc.REASON:
            raise ValueError(f"{acc}: sheet reason is {row['reason_code']!r}, not the "
                             "unmappable-columns pending state; refusing to overwrite")
        if rev["decision"] == "ineligible":
            row.update(verdict="exclude", reason_code=rev["reason_code"], reason_source=SOURCE)
            done["excluded"].append(acc)
            continue
        path = mc.supplementary_file(root, acc)
        if path is None:
            raise FileNotFoundError(f"{acc}: no count file in the cache")
        table = build_table(mc.read_table(path), rev)
        n_samples = table.shape[1] - 1
        if n_samples != int(rev["n_case"]) + int(rev["n_control"]):
            raise ValueError(f"{acc}: {n_samples} columns built, review states "
                             f"{rev['n_case']} + {rev['n_control']}")
        row.update(verdict="pending", reason_code=PENDING, reason_source=SOURCE,
                   annotation_status="annotated", n_case=rev["n_case"],
                   n_control=rev["n_control"])
        done["mapped"].append(acc)
        if write:
            groups = root / rv.CACHE / f"{acc}_groups.csv"
            drop = set(json.loads(rev["drop_json"] or "[]"))
            kept = [r for r in rv._read_csv(groups) if r["sample_id"] not in drop]
            _move_aside(root, acc)
            if drop:
                keep = root / ORIGINALS / acc
                keep.mkdir(parents=True, exist_ok=True)
                shutil.copy(groups, keep / groups.name)
                with open(groups, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=["sample_id", "group"])
                    w.writeheader()
                    w.writerows(kept)
            out = root / rv.CACHE / acc / f"{acc}_raw_count_merged.txt.gz"
            table.to_csv(out, sep="\t", index=False, compression="gzip")
    if write and any(done.values()):
        _write_sheet(sheet, rows, fields)
    return done


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=REPO_ROOT)
    p.add_argument("--reviews", type=Path, default=None)
    p.add_argument("--apply", action="store_true")
    a = p.parse_args()
    done: dict[str, Any] = apply(a.root, a.reviews or a.root / mc.REVIEWS, a.apply)
    verb = "" if a.apply else "would be "
    logger.info("%d series %smapped to pending/%s: %s", len(done["mapped"]), verb, PENDING,
                " ".join(done["mapped"]))
    logger.info("%d series %sexcluded: %s", len(done["excluded"]), verb,
                " ".join(done["excluded"]))


if __name__ == "__main__":
    main()
