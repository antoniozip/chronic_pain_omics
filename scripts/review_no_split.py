"""A local interface for scrutinising the series excluded for having no case/control split.

98 GEO series carry `no_case_control_split_in_series_matrix` in
`literature/prisma/transcriptomics_candidates.csv`, every one set by
`auto_annotate_geo` reading sample metadata against a word list. That exclusion
cannot tell "the deposit has no contrast" from "the parser could not read the
labels" (10 of the series audited so far turned out to be the second). This
serves one page, on this machine only, that puts each series' own header text
and sample table in front of a person and lets them say which it is.

    python scripts/review_no_split.py            # then open http://127.0.0.1:8765
    python scripts/review_no_split.py --reviews temp/trial.csv    # dry run

For each series the reviewer picks the metadata field that carries the
contrast, assigns its values (or a regex) to case and control, and records one
of `contrast_exists`, `no_contrast` or `unclear`. Counts come from the server,
so the page and the saved record cannot disagree about what a rule matches.

Decisions go to `literature/prisma/no_split_review.csv`. **This tool never
writes the screening sheet.** A verdict in the sheet that nothing reads is a
no-op that looks like a decision, and one that something reads moves every
pooled number, so applying the review is a separate, deliberate step.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
import tempfile
import threading
import webbrowser
from collections import Counter
from dataclasses import dataclass
from datetime import UTC, datetime
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "pipeline"))

import auto_annotate_geo as aag  # noqa: E402

REASON = "no_case_control_split_in_series_matrix"
SHEET = Path("literature/prisma/transcriptomics_candidates.csv")
AUDIT = Path("literature/prisma/transcriptomics_rescreen_audit.csv")
REVIEWS = Path("literature/prisma/no_split_review.csv")
CACHE = Path("data/raw/geo_cache")
HTML = Path(__file__).with_name("review_no_split.html")

DECISIONS = ("contrast_exists", "no_contrast", "unclear")
REVIEW_COLUMNS = ["accession", "decision", "field", "rule_json", "n_case",
                  "n_control", "note", "reviewed_at"]
MAX_DISTINCT = 80
MAX_BODY = 1_000_000
HEADER_KEYS = ("Series_title", "Series_summary", "Series_overall_design",
               "Series_type", "Series_pubmed_id", "Series_platform_id")

_WRITE_LOCK = threading.Lock()


@dataclass
class Context:
    root: Path
    reviews: Path
    queue: list[dict[str, Any]]

    def known(self, accession: str) -> bool:
        return any(q["accession"] == accession for q in self.queue)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_reviews(path: Path) -> dict[str, dict[str, str]]:
    return {r["accession"]: r for r in _read_csv(path)}


def build_queue(root: Path, reviews: Path) -> list[dict[str, Any]]:
    """The excluded series, each with the sheet row, its re-audit and its review."""
    audit = {r["accession"]: r for r in _read_csv(root / AUDIT)}
    done = load_reviews(reviews)
    queue = []
    for row in _read_csv(root / SHEET):
        if row["reason_code"] != REASON:
            continue
        acc = row["accession"]
        queue.append({
            "accession": acc, "species": row["species"], "year": row["year"],
            "n_samples": row["n_samples"], "title": row["title"], "url": row["url"],
            "audit": audit.get(acc), "audit_status": (audit.get(acc) or {}).get(
                "status", "not_audited"),
            "review": done.get(acc),
        })
    return sorted(queue, key=lambda q: q["accession"])


def _header_text(path: Path) -> dict[str, str]:
    """The series-level text the matrix carries above its sample table."""
    found: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            key, values = aag._parse_row(line)
            if key in HEADER_KEYS:
                found.setdefault(key, []).extend(v for v in values if v)
    return {k: "\n".join(v) for k, v in found.items()}


def _field(fid: str, label: str, values: list[str], key: str | None) -> dict[str, Any]:
    distinct = Counter(values)
    padded = " " + (key or "").lower() + " "
    skipped = [t for t in aag.SKIP_KEY_TOKENS if key and aag._token_present(t, padded)]
    cand = aag._candidate_from_values(
        label, values, key_drives_case=bool(key) and aag._key_is_pain_signal(key or ""))
    return {
        "id": fid, "label": label, "values": values,
        "n_distinct": len(distinct), "constant": len(distinct) < 2,
        "skipped_by": skipped, "parser": {"case": cand.n_case, "control": cand.n_ctrl},
        "distinct": [{"value": v, "n": n, "parser": aag._classify_value(v)}
                     for v, n in distinct.most_common(MAX_DISTINCT)],
    }


@lru_cache(maxsize=64)
def _load(root: str, accession: str) -> dict[str, Any]:
    cache = Path(root) / CACHE
    files = sorted(cache.glob(f"{accession}_series_matrix.txt.gz")) + sorted(
        cache.glob(f"{accession}-GPL*_series_matrix.txt.gz"))
    parsed = [(p, aag.parse_series_matrix(p)) for p in files]
    parsed = [(p, sm) for p, sm in parsed if sm]
    if not parsed:
        raise FileNotFoundError(f"no readable series matrix for {accession}")
    path, sm = max(parsed, key=lambda t: t[1].n)
    fields = [_field("title", "title", sm.titles, None),
              _field("source", "source_name", sm.source_names, None)]
    fields += [_field(f"c{i}", f"char: {key}", vals, key)
               for i, (key, vals) in enumerate(sm.characteristics)]
    return {"accession": accession, "file": path.name, "gsm": sm.gsm_ids,
            "variants": [{"file": p.name, "n": s.n} for p, s in parsed],
            "header": _header_text(path), "fields": fields}


def series_detail(ctx: Context, accession: str) -> dict[str, Any]:
    if not ctx.known(accession):
        raise KeyError(accession)
    return _load(str(ctx.root), accession)


def _matcher(part: dict[str, Any]):
    exact = set(part.get("values") or [])
    rx = re.compile(part["regex"], re.IGNORECASE) if part.get("regex") else None
    return lambda v: v in exact or (rx is not None and rx.search(v) is not None)


def apply_rule(values: list[str], rule: dict[str, Any]) -> dict[str, Any]:
    """Label each sample case/control/conflict/'' from exact values and regexes."""
    is_case, is_ctrl = _matcher(rule.get("case") or {}), _matcher(rule.get("control") or {})
    labels = []
    for v in values:
        c, k = is_case(v), is_ctrl(v)
        labels.append("conflict" if c and k else "case" if c else "control" if k else "")
    n = Counter(labels)
    return {"labels": labels, "n_case": n["case"], "n_control": n["control"],
            "n_conflict": n["conflict"], "n_unassigned": n[""]}


def preview(ctx: Context, accession: str, field_id: str,
            rule: dict[str, Any]) -> dict[str, Any]:
    fields = {f["id"]: f for f in series_detail(ctx, accession)["fields"]}
    if field_id not in fields:
        raise KeyError(field_id)
    try:
        return apply_rule(fields[field_id]["values"], rule)
    except re.error as exc:
        raise ValueError(f"invalid regex: {exc}") from exc


def save_review(ctx: Context, rec: dict[str, Any]) -> dict[str, str]:
    """Validate one decision and upsert it. `contrast_exists` must be reproducible."""
    acc, decision = str(rec.get("accession", "")), rec.get("decision")
    if not ctx.known(acc):
        raise KeyError(acc)
    if decision not in DECISIONS:
        raise ValueError(f"decision must be one of {', '.join(DECISIONS)}")
    row = dict.fromkeys(REVIEW_COLUMNS, "")
    row.update(accession=acc, decision=decision, note=str(rec.get("note", ""))[:2000],
               reviewed_at=datetime.now(UTC).isoformat(timespec="seconds"))
    if decision == "contrast_exists":
        res = preview(ctx, acc, str(rec.get("field", "")), rec.get("rule") or {})
        if res["n_conflict"]:
            raise ValueError(f"{res['n_conflict']} samples match both case and control")
        if min(res["n_case"], res["n_control"]) < aag.MIN_PER_GROUP:
            raise ValueError(f"needs at least {aag.MIN_PER_GROUP} samples per group, "
                             f"got {res['n_case']} case and {res['n_control']} control")
        label = next(f["label"] for f in series_detail(ctx, acc)["fields"]
                     if f["id"] == rec["field"])
        row.update(field=label, rule_json=json.dumps(
            {"field_id": rec["field"], **rec["rule"]}, sort_keys=True),
            n_case=str(res["n_case"]), n_control=str(res["n_control"]))
    with _WRITE_LOCK:
        merged = load_reviews(ctx.reviews)
        merged[acc] = row
        ctx.reviews.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=ctx.reviews.parent, suffix=".tmp")
        with os.fdopen(fd, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=REVIEW_COLUMNS)
            w.writeheader()
            w.writerows(merged[a] for a in sorted(merged))
        os.replace(tmp, ctx.reviews)
    return row


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
            # A page on another site can reach 127.0.0.1; refuse anything that is
            # not addressed to this server by name, or posted from another origin.
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
            if url.path == "/api/series":
                try:
                    acc = parse_qs(url.query).get("acc", [""])[0]
                    return self._json(200, series_detail(ctx, acc))
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
                if self.path == "/api/preview":
                    return self._json(200, preview(
                        ctx, body["accession"], body["field"], body.get("rule") or {}))
                if self.path == "/api/review":
                    return self._json(200, save_review(ctx, body))
                self._json(404, {"error": "not found"})
            except (KeyError, FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
                self._json(400, {"error": str(exc)})

    return Handler


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--root", type=Path, default=REPO_ROOT)
    ap.add_argument("--reviews", type=Path, default=None,
                    help="where decisions are written (default: the tracked review sheet)")
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()
    reviews = args.reviews or args.root / REVIEWS
    ctx = Context(args.root, reviews, build_queue(args.root, reviews))
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(ctx, args.port))
    url = f"http://127.0.0.1:{args.port}"
    print(f"{len(ctx.queue)} series to review; decisions -> {reviews}\n{url}  (Ctrl-C to stop)")
    if not args.no_open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
