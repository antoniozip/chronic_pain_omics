#!/usr/bin/env python
"""Pool permutation runs of the between-model Q null split across processes.

`scripts/permute_between_model_q.R` can be run as several processes that share
one feature subsample (`--feature-seed`) and draw disjoint label shuffles
(distinct `--seed`). Each writes its own table: the observed statistic as
`perm == 0`, then its permutations numbered from 1, each with its seed. This
script joins them into the single table the paper reads, and refuses two ways
a pooled null could be wrong while looking right:

- runs whose observed rows differ were computed on different features or
  inputs, so their permutations are nulls for different statistics;
- a repeated seed is the same shuffle counted twice, which narrows the null
  and moves the p-value without any new information.

`--reference` names a table whose observed row every run must also reproduce,
usually the one being replaced: a mismatch means the inputs moved since it was
written. The pooled permutations are renumbered 1..B in seed order.

Usage:
    python scripts/pool_permutation_runs.py \\
        --parts results/meta/transcriptomics/stratified/permutation_runs/*.csv \\
        --reference results/meta/transcriptomics/stratified/between_model_Q_permutation.csv \\
        --out results/meta/transcriptomics/stratified/between_model_Q_permutation.csv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
from scipy.stats import beta

logger = logging.getLogger(__name__)

COLUMNS = ["perm", "seed", "tested", "hits", "pct"]


def _observed(table: pd.DataFrame, name: str) -> tuple[int, int]:
    obs = table[table["perm"] == 0]
    if len(obs) != 1:
        raise ValueError(f"{name}: expected one observed row (perm == 0), found {len(obs)}")
    return int(obs["tested"].iloc[0]), int(obs["hits"].iloc[0])


def pool_runs(parts: dict[str, pd.DataFrame],
              reference: tuple[int, int] | None = None) -> pd.DataFrame:
    """Join runs into one table: the shared observed row, then every permutation.

    Args:
        parts: Each run's table, keyed by a name used in error messages.
        reference: (tested, hits) that every run's observed row must match.

    Returns:
        The observed row followed by the permutations in seed order, renumbered
        1..B.

    Raises:
        ValueError: A run lacks seeds, runs disagree on the observed row or
            with `reference`, or a seed occurs twice.
    """
    if not parts:
        raise ValueError("no runs to pool")
    # A nullable integer seed: the observed row has none, and a column that is
    # all missing in one run must not decide the pooled column's type.
    # reindex, so a table written before seeds were recorded reads as seedless.
    parts = {name: t.reindex(columns=COLUMNS).astype({"seed": "Int64"})
             for name, t in parts.items()}
    observed = {name: _observed(t, name) for name, t in parts.items()}
    expected = reference if reference is not None else next(iter(observed.values()))
    wrong = {name: o for name, o in observed.items() if o != expected}
    if wrong:
        raise ValueError(f"observed row (tested, hits) differs from {expected}: {wrong}; "
                         "these runs did not test the same statistic")

    perms = pd.concat([t[t["perm"] > 0].assign(run=name) for name, t in parts.items()],
                      ignore_index=True)
    if perms["seed"].isna().any():
        missing = sorted(perms.loc[perms["seed"].isna(), "run"].unique())
        raise ValueError(f"permutations without a seed in {missing}; "
                         "a repeated shuffle could not be detected")
    repeated = perms[perms["seed"].duplicated(keep=False)]
    if not repeated.empty:
        raise ValueError(f"seed(s) {sorted(repeated['seed'].astype(int).unique())} "
                         f"occur in more than one permutation: {sorted(repeated['run'].unique())}")

    perms = perms.sort_values("seed").reset_index(drop=True)
    perms["perm"] = range(1, len(perms) + 1)
    first = next(iter(parts.values()))
    pooled = pd.concat([first[first["perm"] == 0], perms[COLUMNS]], ignore_index=True)
    return pooled.astype({"perm": int, "tested": int, "hits": int})


def permutation_p(pooled: pd.DataFrame) -> tuple[int, int, float, float, float]:
    """Exceedances b of B, p = (b + 1)/(B + 1), and the Clopper-Pearson 95% interval.

    The interval is for the tail probability b/B estimates: the uncertainty
    that comes from drawing B permutations rather than all of them.
    """
    null = pooled.loc[pooled["perm"] > 0, "pct"]
    observed = float(pooled.loc[pooled["perm"] == 0, "pct"].iloc[0])
    b, n = int((null >= observed).sum()), len(null)
    lo = float(beta.ppf(0.025, b, n - b + 1)) if b > 0 else 0.0
    hi = float(beta.ppf(0.975, b + 1, n - b)) if b < n else 1.0
    return b, n, (b + 1) / (n + 1), lo, hi


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parts", type=Path, nargs="+", required=True)
    p.add_argument("--reference", type=Path, default=None,
                   help="table whose observed row every run must reproduce")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    reference = (_observed(pd.read_csv(args.reference), str(args.reference))
                 if args.reference else None)
    parts = {str(path): pd.read_csv(path) for path in sorted(args.parts)}
    pooled = pool_runs(parts, reference)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pooled.to_csv(args.out, index=False)

    b, n, pval, lo, hi = permutation_p(pooled)
    null = pooled.loc[pooled["perm"] > 0, "pct"]
    logger.info("%d runs, %d permutations -> %s", len(parts), n, args.out)
    logger.info("observed %.2f%%; null mean %.2f%%, 95th percentile %.2f%%",
                float(pooled["pct"].iloc[0]), null.mean(), null.quantile(0.95))
    logger.info("exceeded by %d of %d: p = %.4f (Monte Carlo 95%% CI for the tail "
                "probability %.3f-%.3f)", b, n, pval, lo, hi)


if __name__ == "__main__":
    main()
