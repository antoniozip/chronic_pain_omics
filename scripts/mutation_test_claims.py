#!/usr/bin/env python
"""Mutation-test the manuscript claim registry: does every claim actually fail?

A claim that cannot fail is worse than no claim, because it reports green. The
first version of the gene-claim tests was exactly that: `pytest.approx` keeps a
default absolute tolerance of 1e-12 and passes when *either* tolerance is met,
so 9 of 16 p-value claims compared equal to anything below 1e-12. Four claims
had been spot-checked by hand and looked fine; the rest were vacuous.

This walks the whole registry. For each claim it perturbs the number that claim
reads out of the document that claim names -- manuscript.tex, or
supplementary_body.tex for the supplementary material -- re-runs that claim
alone, and records
whether the perturbation was caught ("killed") or slipped through ("survived").
Survivors are the output that matters.

Every document is restored after every mutation and verified by hash at the
end, so a crash cannot leave one corrupted.

Usage:
    python scripts/mutation_test_claims.py
    python scripts/mutation_test_claims.py --claim pairs_total
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from decimal import Decimal
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tests"))

import test_manuscript_claims as tmc  # noqa: E402

DOCUMENTS = tmc.DOCUMENTS

# Claims that cannot be killed for a stated reason, rather than a defect.
KNOWN_UNKILLABLE = {
    "table_gm35167": (
        "padj underflows to exactly 0, which satisfies any bound; the claim "
        "confirms the direction of the inequality, not the exponent"
    ),
}

# Must cover every number word the manuscript can state, or the claim reading
# it is reported "not perturbable" and checks nothing -- the same hole that let
# "exactly one is---Cfap68" survive its own mutation. It stopped at ten, so
# "Nineteen features are significant" went unverified on 2026-08-30. Kept in
# step with _WORDS in tests/test_manuscript_claims.py.
_WORD_SWAP = {"one": "two", "two": "three", "three": "two", "four": "five",
              "five": "four", "six": "seven", "seven": "six",
              "eight": "nine", "nine": "eight", "ten": "nine",
              "eleven": "twelve", "twelve": "eleven",
              "thirteen": "fourteen", "fourteen": "thirteen",
              "fifteen": "sixteen", "sixteen": "fifteen",
              "seventeen": "eighteen", "eighteen": "seventeen",
              "nineteen": "twenty", "twenty": "nineteen"}


def _perturb(token: str) -> str | None:
    """Change a number as written, minimally but enough to exceed tolerances."""
    t = token.strip()
    if t.lower() in _WORD_SWAP:
        return _WORD_SWAP[t.lower()]
    # A negative is written `$-$0.002`, LaTeX's math minus. Perturb the
    # magnitude and put the sign back: without this the claim reading it is
    # reported "not perturbable" and checks nothing, which is exactly the hole
    # rat's Pearson r sat in while its sign was wrong.
    if t.startswith("$-$"):
        inner = _perturb(t[3:])
        return f"$-${inner}" if inner is not None else None
    if re.fullmatch(r"[\d,]+", t):
        had_commas = "," in t
        value = int(t.replace(",", "")) + 1
        return f"{value:,}" if had_commas else str(value)
    if re.fullmatch(r"-?\d+", t):                       # signed exponent
        return str(int(t) + 1)
    if re.fullmatch(r"\d*\.\d+", t):
        step = Decimal(1).scaleb(-len(t.split(".")[1]))
        return str(Decimal(t) + step)
    return None


def _mutated_text(text: str, claim) -> tuple[str, str, str] | None:
    """Return (mutated text, original token, mutated token) for one claim."""
    match = re.search(claim.pattern, text)
    if match is None:
        return None
    groups = match.groups()
    # Mutate the first group that carries a number we know how to change.
    for idx, token in enumerate(groups, start=1):
        new = _perturb(token)
        if new is None:
            continue
        start, end = match.span(idx)
        return text[:start] + new + text[end:], token, new
    return None


def _run_claim(claim, is_gene: bool) -> str:
    """Run one claim. Returns "passed", "failed" or "skipped".

    A skip means the claim's results file is absent, which is neither a kill
    nor a survival: nothing was checked. Counting it as either would make a
    data-less checkout look verified, or look broken.
    """
    tmc._document.cache_clear()
    try:
        if is_gene:
            tmc.test_named_gene_statistic_matches_results(claim)
        else:
            tmc.test_manuscript_claim_matches_results(claim)
    except AssertionError:
        return "failed"
    except BaseException as exc:
        # pytest's Skipped derives from BaseException, not Exception, so a bare
        # `except Exception` lets it escape and crashes the run on any checkout
        # without results/ -- which is every fresh clone.
        if exc.__class__.__name__ == "Skipped":
            return "skipped"
        raise
    return "passed"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claim", help="only this claim name")
    parser.add_argument(
        "--require-results", action="store_true",
        help="fail if any claim skipped for want of its results file; use in CI "
             "jobs that are supposed to have the data",
    )
    args = parser.parse_args()

    missing = [name for name, path in DOCUMENTS.items() if not path.exists()]
    if missing:
        # Both documents are tracked, so this is a broken checkout rather than
        # the absent-results case the claims themselves skip over.
        print(f"missing document(s): {', '.join(missing)}")
        return 2

    originals = {name: path.read_text() for name, path in DOCUMENTS.items()}
    hashes = {name: hashlib.sha256(text.encode()).hexdigest()
              for name, text in originals.items()}

    # Gene claims are manuscript-only by construction: no p-value is quoted
    # outside it, so GeneClaim carries no source field to read here.
    registry = ([(c, False, c.source) for c in tmc.CLAIMS]
                + [(c, True, "manuscript") for c in tmc.GENE_CLAIMS])
    if args.claim:
        registry = [entry for entry in registry if entry[0].name == args.claim]
        if not registry:
            print(f"no claim named {args.claim!r}")
            return 2

    killed, survived, unmutatable, skipped = [], [], [], []
    try:
        for claim, is_gene, source in registry:
            document = DOCUMENTS[source]
            mutation = _mutated_text(originals[source], claim)
            if mutation is None:
                unmutatable.append(claim.name)
                print(f"  ?  {claim.name}: could not perturb (pattern or token)")
                continue

            text, before, after = mutation
            document.write_text(text)
            outcome = _run_claim(claim, is_gene)
            document.write_text(originals[source])

            if outcome == "skipped":
                skipped.append(claim.name)
                print(f"  -- {claim.name}: skipped   (results not generated)")
            elif outcome == "passed":
                survived.append(claim.name)
                print(f"  !! {claim.name}: SURVIVED  ({before} -> {after})")
            else:
                killed.append(claim.name)
                print(f"  ok {claim.name}: killed    ({before} -> {after})")
    finally:
        for name, path in DOCUMENTS.items():
            path.write_text(originals[name])
        tmc._document.cache_clear()

    for name, path in DOCUMENTS.items():
        restored = hashlib.sha256(path.read_text().encode()).hexdigest()
        assert restored == hashes[name], f"{path.name} was not restored cleanly"

    real_survivors = [n for n in survived if n not in KNOWN_UNKILLABLE]
    print(
        f"\n{len(killed)} killed, {len(survived)} survived "
        f"({len(survived) - len(real_survivors)} expected), "
        f"{len(skipped)} skipped, {len(unmutatable)} not perturbable"
    )
    for name in survived:
        if name in KNOWN_UNKILLABLE:
            print(f"  expected survivor {name}: {KNOWN_UNKILLABLE[name]}")

    if skipped:
        message = (
            f"{len(skipped)} claim(s) were not checked because their results "
            "files are absent. Run the pipeline for real coverage."
        )
        if args.require_results:
            print(f"\nERROR: {message}")
            return 1
        print(f"\nWARNING: {message}")

    if real_survivors or unmutatable:
        print("\nclaims that do not detect a wrong number:")
        for name in real_survivors + unmutatable:
            print(f"  - {name}")
        return 1

    if skipped:
        # Do not report a clean bill of health over a handful of claims when
        # most of the registry was never exercised.
        print(f"the {len(killed)} claim(s) that ran all detect a wrong number; "
              f"{len(skipped)} were not checked.")
    else:
        print("every claim detects a wrong number.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
