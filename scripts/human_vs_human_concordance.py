#!/usr/bin/env python3
"""Positive controls for the cross-species concordance analysis.

The paper's central cross-species result is a null: rodent-to-human directional
concordance sits at chance. A null is only interpretable if the machinery is
shown capable of detecting agreement that is really there, and nothing in the
analysis demonstrated that.

This measures concordance between independent halves of the *human* corpus,
pooled through the same ``06g_species_meta.R`` code path as the published
pools, under three designs of decreasing difficulty:

    random      five random halves of all 23 human studies (mixed compartment
                on both sides)
    compartment the 13 endometrial studies against the other 10 (same species,
                mismatched compartment)
    replicate   the endometriosis stratum split against itself, three seeds
                (same condition, same compartment, same species -- the cleanest
                replicate pair the corpus can offer)

Halves are joined directly on gene symbol, with none of the ortholog-mapping
loss the cross-species join carries, so every figure here is an *upper bound*
on what the machinery can detect.

Concordance is defined exactly as in
``src/cp_multiomics/ortholog/concordance.py``: the share of jointly measured
features whose pooled effects share a sign.

Result (2026-09-03): every design returns 50-53%, and none yields a single
feature significant in both halves. The cross-species figures (49-51%) lie
inside the range of same-condition replicates, so the cross-species comparison
cannot separate a fact about translation from the reproducibility floor of the
measurement. See plan/2026-09-03-small-study-effects-rebuttal.md.

The halves themselves are built by --build-halves, which derives every split
from the *current* human pool and writes one exclusion list per half. Until
2026-09-12 they came from ad-hoc shell scripts under temp/hvh, and the lists
were fixed at the 23-study pool they were written for. Against the 34-study
pool those lists leave 11 studies in *both* halves of every split, which
destroys the independence the control exists to demonstrate -- and the script
would have reported a concordance figure anyway, since nothing checked that
the halves were disjoint. They are now checked.

Usage:
    python scripts/human_vs_human_concordance.py --build-halves
    python scripts/human_vs_human_concordance.py            # compare only
"""

from __future__ import annotations

import argparse
import logging
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

#: (design, label, directory A, directory B). Directories are 06g --out-dir
#: parents produced by temp/hvh/run*.sh.
COMPARISONS: list[tuple[str, str, str, str]] = (
    [("random", f"seed {s}", f"s{s}A", f"s{s}B") for s in range(1, 6)]
    + [("compartment", "endometrial vs rest", "s0E", "s0N")]
    + [("replicate", f"seed {s}", f"e{s}A", f"e{s}B") for s in (7, 8, 9)]
)


COMPARTMENTS = REPO_ROOT / "conf" / "analysis" / "human_tissue_compartments.csv"
POOLED_HUMAN = (REPO_ROOT / "results" / "meta" / "transcriptomics"
                / "by_species" / "Homo_sapiens_pooled.csv")


def human_units() -> list[str]:
    """Study units in the unrestricted human pool, as 06g wrote it."""
    if not POOLED_HUMAN.exists():
        raise SystemExit(f"run 06g first: {POOLED_HUMAN} is missing")
    d = pd.read_csv(POOLED_HUMAN, usecols=["study_ids"])
    return sorted(d["study_ids"].astype(str).str.split(";")
                  .explode().str.strip().unique())


def endometrial_units(units: list[str]) -> list[str]:
    d = pd.read_csv(COMPARTMENTS, dtype=str).set_index("study_id")["compartment"]
    missing = [u for u in units if u.split("_")[0] not in d.index]
    if missing:
        raise SystemExit(
            f"no compartment for: {', '.join(missing)} -- annotate them in "
            f"{COMPARTMENTS.relative_to(REPO_ROOT)}")
    return [u for u in units if d.loc[u.split("_")[0]] == "endometrial"]


def planned_halves() -> dict[str, list[str]]:
    """Each half's kept studies, derived from the current pool.

    Random splits shuffle all human units; the compartment split is the
    endometrial studies against the rest; the replicate split halves the
    endometrial set, which is the cleanest same-condition pair the corpus can
    offer.
    """
    units = human_units()
    endo = endometrial_units(units)
    rest = [u for u in units if u not in set(endo)]
    keep: dict[str, list[str]] = {}
    for seed in range(1, 6):
        shuffled = list(units)
        random.Random(seed).shuffle(shuffled)
        mid = len(shuffled) // 2
        keep[f"s{seed}A"] = sorted(shuffled[:mid])
        keep[f"s{seed}B"] = sorted(shuffled[mid:])
    keep["s0E"], keep["s0N"] = sorted(endo), sorted(rest)
    for seed in (7, 8, 9):
        shuffled = list(endo)
        random.Random(seed).shuffle(shuffled)
        mid = len(shuffled) // 2
        keep[f"e{seed}A"] = sorted(shuffled[:mid])
        keep[f"e{seed}B"] = sorted(shuffled[mid:])
    return keep


def build_halves(work_dir: Path) -> None:
    """Pool each half through 06g, the same code path the published pools use."""
    units = human_units()
    keep = planned_halves()
    for design, _label, a, b in COMPARISONS:
        overlap = set(keep[a]) & set(keep[b])
        if overlap:
            raise SystemExit(
                f"{design} halves {a}/{b} share {len(overlap)} studies: "
                f"{', '.join(sorted(overlap))}. A positive control whose halves "
                "overlap measures nothing.")
    work_dir.mkdir(parents=True, exist_ok=True)
    for name, kept in keep.items():
        out = work_dir / name / "by_species"
        if (out / "Homo_sapiens_pooled.csv").exists():
            logger.info("%s already built; skipping", name)
            continue
        excluded = sorted(set(units) - set(kept))
        ex = work_dir / f"exclude_{name}.csv"
        ex.write_text("study_id\n" + "\n".join(excluded) + "\n")
        out.mkdir(parents=True, exist_ok=True)
        cmd = ["Rscript", "pipeline/06g_species_meta.R",
               "--modality", "transcriptomics",
               "--exclude-studies", str(ex), "--out-dir", str(out)]
        logger.info("%s: %d studies kept", name, len(kept))
        r = subprocess.run(cmd, cwd=REPO_ROOT,
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            raise SystemExit(f"06g failed for {name}")


def half_pool(work_dir: Path, name: str) -> pd.DataFrame:
    path = work_dir / name / "by_species" / "Homo_sapiens_pooled.csv"
    if not path.exists():
        raise SystemExit(
            f"missing half-pool: {path}\nBuild them first: "
            "python scripts/human_vs_human_concordance.py --build-halves")
    return pd.read_csv(path, usecols=["feature_id", "k", "yi", "se", "padj"])


def compare(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Concordance between two independent pools of the same species."""
    m = a.merge(b, on="feature_id", suffixes=("_a", "_b")).dropna(
        subset=["yi_a", "yi_b"])
    if len(m) < 3:
        raise SystemExit("fewer than three jointly measured features")

    conc = np.sign(m["yi_a"]) == np.sign(m["yi_b"])
    n, c = len(m), int(conc.sum())
    rho, rho_p = stats.spearmanr(m["yi_a"], m["yi_b"])

    kmin = m[["k_a", "k_b"]].min(axis=1)
    deep = m[kmin >= 3]
    deep_pct = (100.0 * float((np.sign(deep["yi_a"]) == np.sign(deep["yi_b"])).mean())
                if len(deep) else float("nan"))

    # The diagnostic that decides whether the measurement is sound: on a sound
    # one the best-estimated features are the most reproducible.
    semax = m[["se_a", "se_b"]].max(axis=1)
    precise = m[semax <= semax.quantile(0.10)]
    precise_pct = (100.0 * float(
        (np.sign(precise["yi_a"]) == np.sign(precise["yi_b"])).mean())
        if len(precise) else float("nan"))

    return {
        "n_features": n,
        "n_concordant": c,
        "pct_concordant": 100.0 * c / n,
        "binomial_p": float(stats.binomtest(c, n, 0.5).pvalue),
        "spearman_r": float(rho),
        "spearman_p": float(rho_p),
        "pct_concordant_k3": deep_pct,
        "pct_concordant_precise10": precise_pct,
        "n_both_significant": int(
            ((m["padj_a"] < 0.05) & (m["padj_b"] < 0.05)).sum()),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", type=Path, default=REPO_ROOT / "temp" / "hvh")
    ap.add_argument("--build-halves", action="store_true",
                    help="pool each half through 06g before comparing")
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "results" / "cross_species"
                    / "transcriptomics" / "human_vs_human.csv")
    args = ap.parse_args()

    if args.build_halves:
        build_halves(args.work_dir)

    rows = []
    for design, label, a, b in COMPARISONS:
        res = {"design": design, "label": label}
        res.update(compare(half_pool(args.work_dir, a),
                           half_pool(args.work_dir, b)))
        rows.append(res)
        logger.info(
            "%-11s %-20s n=%6d  %.2f%%  (k>=3 %.2f%%, precise %.2f%%)  "
            "rho=%+.3f  both-sig=%d",
            design, label, res["n_features"], res["pct_concordant"],
            res["pct_concordant_k3"], res["pct_concordant_precise10"],
            res["spearman_r"], res["n_both_significant"],
        )

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    print()
    for design, g in df.groupby("design", sort=False):
        logger.info("%-11s concordance %.2f%% (%.2f--%.2f), both-significant %d",
                    design, g["pct_concordant"].mean(), g["pct_concordant"].min(),
                    g["pct_concordant"].max(), int(g["n_both_significant"].sum()))
    logger.info("-> %s", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
