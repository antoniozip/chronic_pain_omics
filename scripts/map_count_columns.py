"""A local interface for tying a count file's columns to the samples GEO describes.

Eleven series sit in `literature/prisma/transcriptomics_candidates.csv` as
`pending / count_columns_not_mappable_to_samples`: their design was cleared by
hand (`no_split_review.csv`, one `{acc}_groups.csv` each), but step 05 cannot say
which column of the deposited count file is which sample. The columns carry the
authors' own labels (`NE_1`, `CPIP1 FPKM`, `N381_b_1`), the series matrix carries
GSM ids and titles, and nothing joins them. This serves one page, on this machine
only, that proposes the join, shows each pairing beside the sample's title and its
case/control label, and records what a person confirms.

    python scripts/map_count_columns.py          # then open http://127.0.0.1:8766
    python scripts/map_count_columns.py --reviews temp/trial.csv    # dry run

A pairing is proposed only when it is unambiguous (label equal, contained, or a
subsequence of the sample's title or library name, and claimed by one column
alone); columns left over are paired by position and marked as such. Three
series were merged by step 05d from a groups file, which names its columns by
group and index; the 05d log lists the GSM behind each, and that is read back.

Decisions go to `literature/prisma/count_column_mapping.csv`. **This tool writes
neither the screening sheet nor the cache**; `apply_count_mapping.py` does, as a
separate step.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import threading
import webbrowser
from dataclasses import dataclass
from datetime import UTC, datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import review_no_split as rv  # noqa: E402

REASON = "count_columns_not_mappable_to_samples"
REVIEWS = Path("literature/prisma/count_column_mapping.csv")
LOG_05D = Path("temp/run_05d.log")
HTML = Path(__file__).with_name("map_count_columns.html")

DECISIONS = ("mapping_ok", "ineligible", "unclear")
INELIGIBLE_REASONS = ("in_vitro_stimulus_not_patient_contrast",
                      "in_vitro_hormone_treatment_not_patient_contrast",
                      "cultured_cells_across_passage_and_cell_type",
                      "no_pain_versus_control_contrast", "no_pain_phenotype",
                      "genotype_contrast_not_pain_versus_control")
SCALES = ("counts", "linear", "log2")
COLUMNS = ["accession", "decision", "reason_code", "gene_col", "scale", "mapping_json",
           "drop_json", "n_case", "n_control", "note", "reviewed_at"]
PREVIEW_ROWS = 200
MAX_BODY = 1_000_000
# Words that decorate a column or a title without naming the sample.
NOISE = {"rep", "replicate", "repeat", "sample", "fpkm", "tpm", "norm", "normalized",
         "normalised", "counts", "count", "readcount", "raw", "library", "name"}
_WRITE_LOCK = threading.Lock()


@dataclass
class Context:
    root: Path
    reviews: Path
    queue: list[dict[str, Any]]

    def known(self, accession: str) -> bool:
        return any(q["accession"] == accession for q in self.queue)


# --------------------------------------------------------------------------- matching

def label_key(text: str) -> str:
    """A label reduced to its letters and numbers, without the decorating words."""
    tokens = re.findall(r"[a-z]+|\d+", text.lower())
    return "".join(t for t in tokens if t not in NOISE)


def _digits(key: str) -> list[str]:
    return re.findall(r"\d+", key)


def _subsequence(short: str, long: str) -> bool:
    it = iter(long)
    return all(ch in it for ch in short)


def pair_score(column: str, strings: list[str]) -> tuple[int, str]:
    """How well a column label names a sample: 3 equal, 2 contained, 1 subsequence."""
    col = label_key(column)
    if len(col) < 2:
        return 0, ""
    best = (0, "")
    for s in strings:
        key = label_key(s)
        if not key:
            continue
        if key == col:
            return 3, "equal"
        if col in key or (len(key) >= 3 and key in col):
            best = max(best, (2, "contained"))
        elif (len(col) >= 3 and _subsequence(col, key) and _digits(col) == _digits(key)[-len(
                _digits(col)):] and _digits(col)):
            best = max(best, (1, "subsequence"))
    return best


def propose(columns: list[str], samples: list[dict[str, Any]]) -> dict[str, dict[str, str]]:
    """Pair columns with samples where one pairing is clearly the best, then by position.

    Returns `{column: {"gsm", "how"}}`. A column whose best score is shared by two
    samples, or a sample claimed at that score by two columns, is left to position.
    """
    scores = {c: {s["gsm"]: pair_score(c, s["strings"]) for s in samples} for c in columns}
    out: dict[str, dict[str, str]] = {}
    for c in columns:
        top = max((v[0] for v in scores[c].values()), default=0)
        if top == 0:
            continue
        hits = [g for g, v in scores[c].items() if v[0] == top]
        if len(hits) != 1:
            continue
        gsm = hits[0]
        rivals = [d for d in columns if d != c and scores[d][gsm][0] >= top]
        if not rivals:
            out[c] = {"gsm": gsm, "how": scores[c][gsm][1]}
    left_cols = [c for c in columns if c not in out]
    left_gsm = [s["gsm"] for s in samples if s["gsm"] not in {v["gsm"] for v in out.values()}]
    block = _positional_block(left_cols, len(left_gsm))
    out.update({c: {"gsm": g, "how": "positional"} for c, g in zip(block, left_gsm)})
    return out


def _is_raw_count_label(name: str) -> bool:
    low = name.lower()
    return "count" in low and not ("norm" in low or "average" in low)


def _positional_block(columns: list[str], n: int) -> list[str]:
    """`n` consecutive columns to pair with the leftover samples by position.

    A file that deposits counts, FPKM and TPM side by side has several blocks of
    the right length; the one named as raw read counts is taken, because the
    others cannot be told apart by size. With no such block, only an exact fit
    is paired.
    """
    if n == 0 or len(columns) < n:
        return []
    if len(columns) == n:
        return columns
    windows = [columns[i:i + n] for i in range(len(columns) - n + 1)]
    best = max(windows, key=lambda w: sum(_is_raw_count_label(c) for c in w))
    return best if all(_is_raw_count_label(c) for c in best) else []


def read_05d_log(path: Path, accession: str) -> dict[str, str]:
    """`{column: GSM}` as step 05d logged it for a merged table named by group."""
    if not path.exists():
        return {}
    lines = path.read_text(errors="replace").splitlines()
    try:
        start = lines.index(f"=== {accession} ===")
    except ValueError:
        return {}
    found: dict[str, str] = {}
    for line in lines[start + 1:]:
        if line.startswith("=== "):
            break
        m = re.match(r"\s+(GSM\d+) → (\S+) \(", line)
        if m:
            found[m.group(2)] = m.group(1)
    return found


# --------------------------------------------------------------------------- tables

def supplementary_file(root: Path, accession: str) -> Path | None:
    """The one count file the cache holds for a series, whatever shape it came in."""
    folder = root / rv.CACHE / accession
    files = sorted(p for p in folder.glob("*") if p.is_file()) if folder.is_dir() else []
    merged = [p for p in files if p.name.endswith("_raw_count_merged.txt.gz")]
    return (merged or files or [None])[0]  # type: ignore[list-item]


def read_table(path: Path, nrows: int | None = None) -> pd.DataFrame:
    """A count table with one header row, whichever of three layouts it came in.

    A few deposits put a decorative title in row 1 and the sample names in row 2
    (`Normalized RC (log2)` over `WT_1`, `WT_2`...); those are recognised by a
    second row that is text across the sample columns, and re-headed from it.
    """
    if path.suffix == ".xlsx":
        df = pd.read_excel(path, nrows=nrows)
    else:
        df = pd.read_csv(path, sep="," if ".csv" in path.name else "\t", nrows=nrows)
    if len(df) and df.shape[1] > 2:
        row = df.iloc[0, 1:]
        if row.notna().sum() >= df.shape[1] - 2 and not any(_is_number(v) for v in row):
            df.columns = [df.columns[0], *[str(v) for v in row]]
            df = df.iloc[1:].reset_index(drop=True)
    first = str(df.columns[0])
    if first.startswith("Unnamed"):
        df = df.rename(columns={df.columns[0]: "gene_id"})
    return df


def _is_number(value: Any) -> bool:
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def describe_table(df: pd.DataFrame) -> dict[str, Any]:
    cols = []
    for name in df.columns:
        col = df[name].dropna()
        numeric = len(col) > 0 and col.map(_is_number).mean() > 0.95
        cols.append({"name": str(name), "numeric": bool(numeric),
                     "preview": [str(v)[:24] for v in df[name].head(3)]})
    return {"columns": cols}


def guess_scale(df: pd.DataFrame, columns: list[str]) -> dict[str, Any]:
    """Counts, linear or log2, read off the columns that will be used."""
    nums = df[columns].apply(pd.to_numeric, errors="coerce") if columns else df.iloc[:, :0]
    values = nums.stack() if not nums.empty else pd.Series(dtype=float)
    if not len(values):
        return {"scale_guess": "linear", "min": 0.0, "max": 0.0}
    whole = float(((values % 1) == 0).mean())
    low, top = float(values.min()), float(values.max())
    guess = "counts" if whole > 0.99 else "log2" if top < 50 and low >= 0 else "linear"
    return {"scale_guess": guess, "min": low, "max": top}


def gene_column_guess(columns: list[dict[str, Any]]) -> str:
    names = [c["name"] for c in columns if not c["numeric"]]
    for pat in (r"symbol|gene.?name", r"gene"):
        for n in names:
            if re.search(pat, n, re.I):
                return n
    return names[0] if names else columns[0]["name"]


# --------------------------------------------------------------------------- study

def _groups(root: Path, accession: str) -> dict[str, str]:
    path = root / rv.CACHE / f"{accession}_groups.csv"
    return {r["sample_id"]: r["group"] for r in rv._read_csv(path)}


def _samples(root: Path, accession: str) -> list[dict[str, Any]]:
    detail = rv._load(str(root), accession)
    groups = _groups(root, accession)
    out = []
    for i, gsm in enumerate(detail["gsm"]):
        fields = {f["label"]: f["values"][i] for f in detail["fields"]}
        strings = [fields.get("title", ""), fields.get("source_name", "")]
        out.append({"gsm": gsm, "title": fields.get("title", ""),
                    "group": groups.get(gsm, ""),
                    "details": [f"{k.removeprefix('char: ')}: {v}" for k, v in fields.items()
                                if k.startswith("char:")],
                    "strings": strings})
    return out


def study_detail(ctx: Context, accession: str) -> dict[str, Any]:
    if not ctx.known(accession):
        raise KeyError(accession)
    detail = rv._load(str(ctx.root), accession)
    samples = _samples(ctx.root, accession)
    path = supplementary_file(ctx.root, accession)
    if path is None:
        raise FileNotFoundError(f"no count file cached for {accession}")
    preview = read_table(path, PREVIEW_ROWS)
    table = describe_table(preview)
    numeric = [c["name"] for c in table["columns"] if c["numeric"]]
    logged = read_05d_log(ctx.root / LOG_05D, accession)
    if logged and set(logged) <= set(numeric):
        proposal = {c: {"gsm": g, "how": "05d log"} for c, g in logged.items()}
    else:
        proposal = propose(numeric, samples)
    table.update(guess_scale(preview, list(proposal) or numeric))
    review = rv.load_reviews(ctx.reviews).get(accession)
    return {"accession": accession, "header": detail["header"], "file": path.name,
            "samples": [{k: v for k, v in s.items() if k != "strings"} for s in samples],
            "table": table, "gene_col": gene_column_guess(table["columns"]),
            "proposal": proposal, "review": review}


def build_queue(root: Path, reviews: Path) -> list[dict[str, Any]]:
    done = rv.load_reviews(reviews)
    queue = []
    for row in rv._read_csv(root / rv.SHEET):
        if row["reason_code"] != REASON:
            continue
        acc = row["accession"]
        queue.append({"accession": acc, "species": row["species"], "year": row["year"],
                      "n_samples": row["n_samples"], "title": row["title"],
                      "url": row["url"], "review": done.get(acc)})
    return sorted(queue, key=lambda q: q["accession"])


def validate(ctx: Context, rec: dict[str, Any]) -> dict[str, Any]:
    """Check a decision against the series and return the counts the record will state."""
    acc = rec["accession"]
    if not ctx.known(acc):
        raise KeyError(acc)
    decision = rec.get("decision")
    if decision not in DECISIONS:
        raise ValueError(f"decision must be one of {', '.join(DECISIONS)}")
    if decision == "ineligible":
        if rec.get("reason_code") not in INELIGIBLE_REASONS:
            raise ValueError("an ineligible series needs a reason code")
        if not (rec.get("note") or "").strip():
            raise ValueError("an ineligible series needs a note saying why")
    if decision != "mapping_ok":
        return {"n_case": "", "n_control": ""}
    detail = study_detail(ctx, acc)
    mapping: dict[str, str] = rec.get("mapping") or {}
    drop = set(rec.get("drop") or [])
    columns = {c["name"] for c in detail["table"]["columns"] if c["numeric"]}
    gsm = {s["gsm"]: s["group"] for s in detail["samples"]}
    if rec.get("scale") not in SCALES:
        raise ValueError(f"scale must be one of {', '.join(SCALES)}")
    if rec.get("gene_col") not in {c["name"] for c in detail["table"]["columns"]}:
        raise ValueError("the gene column is not in the file")
    if rec["gene_col"] in mapping:
        raise ValueError("the gene column cannot also be a sample column")
    bad = sorted(c for c in mapping if c not in columns)
    if bad:
        raise ValueError(f"not numeric columns of the file: {bad[:3]}")
    unknown = sorted(g for g in [*mapping.values(), *drop] if g not in gsm)
    if unknown:
        raise ValueError(f"not samples of this series: {unknown[:3]}")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("two columns are paired with the same sample")
    wanted = {g for g, grp in gsm.items() if grp in ("case", "control")} - drop
    if wanted - set(mapping.values()):
        raise ValueError(f"{len(wanted - set(mapping.values()))} case/control samples have no "
                         "column; pair them or drop them from the groups")
    arms = [sum(gsm[g] == arm for g in mapping.values() if g not in drop)
            for arm in ("case", "control")]
    if min(arms) < 2:
        raise ValueError(f"at least 2 per arm are needed, got {arms[0]} case / {arms[1]} control")
    return {"n_case": str(arms[0]), "n_control": str(arms[1])}


def save_review(ctx: Context, rec: dict[str, Any]) -> dict[str, str]:
    counts = validate(ctx, rec)
    row = {"accession": rec["accession"], "decision": rec["decision"],
           "reason_code": rec.get("reason_code") or "", "gene_col": rec.get("gene_col") or "",
           "scale": rec.get("scale") or "",
           "mapping_json": json.dumps(rec.get("mapping") or {}, sort_keys=True),
           "drop_json": json.dumps(sorted(rec.get("drop") or [])),
           "note": (rec.get("note") or "").strip(),
           "reviewed_at": datetime.now(UTC).isoformat(timespec="seconds"), **counts}
    with _WRITE_LOCK:
        rows = rv.load_reviews(ctx.reviews)
        rows[row["accession"]] = row
        ctx.reviews.parent.mkdir(parents=True, exist_ok=True)
        tmp = ctx.reviews.with_suffix(".tmp")
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(rows[k] for k in sorted(rows))
        tmp.replace(ctx.reviews)
    return row


# --------------------------------------------------------------------------- server

def make_handler(ctx: Context, port: int) -> type[BaseHTTPRequestHandler]:
    hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args: Any) -> None:
            pass

        def _send(self, code: int, body: bytes, ctype: str) -> None:
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, code: int, obj: Any) -> None:
            self._send(code, json.dumps(obj).encode(), "application/json")

        def _local(self) -> bool:
            origin = self.headers.get("Origin")
            return self.headers.get("Host") in hosts and (
                origin is None or origin.removeprefix("http://") in hosts)

        def do_GET(self) -> None:
            if not self._local():
                return self._json(403, {"error": "forbidden"})
            url = urlparse(self.path)
            if url.path == "/":
                return self._send(200, HTML.read_bytes(), "text/html; charset=utf-8")
            if url.path == "/api/queue":
                return self._json(200, build_queue(ctx.root, ctx.reviews))
            if url.path == "/api/study":
                try:
                    return self._json(200, study_detail(
                        ctx, parse_qs(url.query).get("acc", [""])[0]))
                except (KeyError, FileNotFoundError) as exc:
                    return self._json(404, {"error": str(exc)})
            self._json(404, {"error": "not found"})

        def do_POST(self) -> None:
            if not self._local() or self.headers.get("Content-Type") != "application/json":
                return self._json(403, {"error": "forbidden"})
            size = int(self.headers.get("Content-Length") or 0)
            if size > MAX_BODY:
                return self._json(413, {"error": "too large"})
            try:
                body = json.loads(self.rfile.read(size) or b"{}")
                if self.path == "/api/review":
                    return self._json(200, save_review(ctx, body))
                self._json(404, {"error": "not found"})
            except (KeyError, FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
                self._json(400, {"error": str(exc)})

    return Handler


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--port", type=int, default=8766)
    ap.add_argument("--root", type=Path, default=REPO_ROOT)
    ap.add_argument("--reviews", type=Path, default=None)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()
    reviews = args.reviews or args.root / REVIEWS
    ctx = Context(args.root, reviews, build_queue(args.root, reviews))
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(ctx, args.port))
    url = f"http://127.0.0.1:{args.port}"
    print(f"{len(ctx.queue)} series to map; decisions -> {reviews}\n{url}  (Ctrl-C to stop)")
    if not args.no_open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
