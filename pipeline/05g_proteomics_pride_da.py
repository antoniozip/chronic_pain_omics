#!/usr/bin/env python
"""Per-study differential abundance for the PRIDE proteomics studies (Phase 3).

Emits the same schema `06_meta_analysis.R` consumes for every other modality,
keyed on **gene symbol** so the studies can pool with one another:

    feature_id, effect_size, se, pval, control_mean, treatment_mean,
    n_control, n_treatment, study_id, region, condition, cohort, padj

`effect_size` is Hedges' g, matching the existing PXD013362 arm.

Contrasts (settled by reading each table's design, see
plan/2026-08-14-proteomics-expansion.md):

  PXD054342  cci vs sham. The `ex` prior-exercise arm is excluded rather than
             folded into either side: in the case group it would measure
             exercise, in the control group it would measure protection.
  PXD055816  SNI vs Sham, naive dropped, pooling sex/age/timepoint into one
             estimate per feature as every other study here contributes one.
             Pooling D7 and D14 averages an acute and a persistent phase.

Usage:
    python pipeline/05g_proteomics_pride_da.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
TABLES = REPO_ROOT / "data" / "raw" / "proteomics" / "tables"
OUT_ROOT = REPO_ROOT / "results" / "per_study" / "proteomics"

MIN_PER_GROUP = 3          # usable observations required in each arm
SCHEMA = ["feature_id", "effect_size", "se", "pval", "control_mean",
          "treatment_mean", "n_control", "n_treatment", "study_id",
          "region", "condition", "cohort", "padj"]


def hedges_g(case: np.ndarray, ctrl: np.ndarray) -> tuple[float, float]:
    """Hedges' g and its standard error for two independent samples."""
    n1, n2 = len(case), len(ctrl)
    df = n1 + n2 - 2
    if df <= 0:
        return np.nan, np.nan
    v1, v2 = case.var(ddof=1), ctrl.var(ddof=1)
    sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    if not np.isfinite(sp2) or sp2 <= 0:
        return np.nan, np.nan
    d = (case.mean() - ctrl.mean()) / np.sqrt(sp2)
    j = 1.0 - 3.0 / (4.0 * df - 1.0)          # small-sample correction
    g = j * d
    var_g = (j ** 2) * ((n1 + n2) / (n1 * n2) + d ** 2 / (2.0 * df))
    return g, float(np.sqrt(var_g))


def differential(expr: pd.DataFrame, case_cols: list[str], ctrl_cols: list[str],
                 study_id: str, condition: str, region: str) -> pd.DataFrame:
    rows = []
    case_mat, ctrl_mat = expr[case_cols].to_numpy(), expr[ctrl_cols].to_numpy()
    for i, feature in enumerate(expr.index):
        case = case_mat[i][np.isfinite(case_mat[i])]
        ctrl = ctrl_mat[i][np.isfinite(ctrl_mat[i])]
        if len(case) < MIN_PER_GROUP or len(ctrl) < MIN_PER_GROUP:
            continue
        g, se = hedges_g(case, ctrl)
        if not np.isfinite(g) or not np.isfinite(se) or se <= 0:
            continue
        p = stats.ttest_ind(case, ctrl, equal_var=False).pvalue
        rows.append({
            "feature_id": feature, "effect_size": g, "se": se, "pval": p,
            "control_mean": ctrl.mean(), "treatment_mean": case.mean(),
            "n_control": len(ctrl), "n_treatment": len(case),
            "study_id": study_id, "region": region, "condition": condition,
            "cohort": "", "padj": np.nan,
        })
    df = pd.DataFrame(rows, columns=SCHEMA)
    if not df.empty:
        from statsmodels.stats.multitest import multipletests
        df["padj"] = multipletests(df["pval"], method="fdr_bh")[1]
    return df


def to_log2(frame: pd.DataFrame) -> pd.DataFrame:
    """MaxQuant writes absent measurements as 0; they are missing, not zero."""
    out = frame.apply(pd.to_numeric, errors="coerce")
    return np.log2(out.where(out > 0))


def collapse_to_gene(frame: pd.DataFrame, genes: pd.Series) -> pd.DataFrame:
    """One row per gene symbol, taking the most intense entry per gene.

    Protein groups map many-to-one onto genes; summing log2 intensities would
    be meaningless and averaging them would mix isoform-level noise, so the
    best-measured group represents the gene.
    """
    key = genes.fillna("").str.split(";").str[0].str.strip()
    frame = frame[key != ""]
    key = key[key != ""]
    order = frame.mean(axis=1).fillna(-np.inf)
    keep = order.groupby(key).idxmax()
    out = frame.loc[keep]
    out.index = keep.index
    out.index.name = "feature_id"
    return out


def pxd054342() -> pd.DataFrame:
    path = TABLES / "PXD054342_proteinGroups.txt"
    df = pd.read_csv(path, sep="\t", low_memory=False)
    for flag in ("Reverse", "Potential contaminant", "Only identified by site"):
        if flag in df.columns:
            df = df[df[flag].isna() | (df[flag].astype(str).str.strip() != "+")]
    cols = [c for c in df.columns if c.startswith("Intensity ")]
    expr = collapse_to_gene(to_log2(df[cols]), df["Gene names"])
    case = [c for c in expr.columns if re.match(r"Intensity cci-", c)]
    ctrl = [c for c in expr.columns if re.match(r"Intensity sham-", c)]
    print(f"  PXD054342: {expr.shape[0]} genes | cci={len(case)} sham={len(ctrl)} "
          f"(excluded {len([c for c in cols if 'ex-' in c])} prior-exercise samples)")
    return differential(expr, case, ctrl, "PXD054342", "CCI", "spinal_dorsal_horn")


def pxd055816() -> pd.DataFrame:
    path = TABLES / "PXD055816_GeneGroup_Quant.csv"
    df = pd.read_csv(path, low_memory=False)
    gene_col = df.columns[0]
    samples = [c for c in df.columns[1:]]
    expr = to_log2(df[samples])
    expr.index = df[gene_col].astype(str).str.strip()
    expr = expr[expr.index != ""]
    expr = expr.groupby(level=0).max()
    expr.index.name = "feature_id"
    case = [c for c in samples if c.startswith("SNI_")]
    ctrl = [c for c in samples if c.startswith("Sham_")]
    naive = [c for c in samples if c.startswith("NV_")]
    print(f"  PXD055816: {expr.shape[0]} genes | SNI={len(case)} Sham={len(ctrl)} "
          f"(excluded {len(naive)} naive)")
    return differential(expr, case, ctrl, "PXD055816", "SNI", "DRG")


def main() -> int:
    print("Per-study proteomics DA (Hedges' g, gene-symbol keyed)")
    for accession, build in (("PXD054342", pxd054342), ("PXD055816", pxd055816)):
        table = build()
        if table.empty:
            print(f"  {accession}: no feature met the {MIN_PER_GROUP}-per-group minimum")
            continue
        out_dir = OUT_ROOT / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f"{accession}_effects.csv"
        table.to_csv(out, index=False)
        sig = int((table["padj"] < 0.05).sum())
        print(f"  -> {out.relative_to(REPO_ROOT)}  {len(table)} features, "
              f"{sig} at padj<0.05")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
