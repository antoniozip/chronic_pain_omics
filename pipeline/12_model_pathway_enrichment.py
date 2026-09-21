"""Per-pain-model pathway enrichment (Step 4.2).

Runs Enrichr ORA for each pain model subgroup (CCI, SNL) using genes
significant within that model's stratified meta-analysis (06b outputs).
SNI (1 sig gene) and CFA (8 sig genes) are reported descriptively.

Inputs:
  results/meta/transcriptomics/stratified/{model}_pooled.csv
  data/interim/ortholog_cache.json

Outputs:
  results/pathway_enrichment/models/
    {model}_{library}.csv
    {model}_{library}_top20.pdf
    model_comparison.csv   — top enriched terms per model side by side
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

REPO_ROOT  = Path(__file__).resolve().parent.parent
STRAT_DIR  = REPO_ROOT / "results" / "meta" / "transcriptomics" / "stratified"
CACHE_JSON = REPO_ROOT / "data" / "interim" / "ortholog_cache.json"
OUT_DIR    = REPO_ROOT / "results" / "pathway_enrichment" / "models"

PADJ_THRESH = 0.05
K_MIN       = 2          # require k >= 2 cross-study replication
MIN_LIST    = 10         # minimum genes for Enrichr ORA

ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
LIBRARIES   = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "WikiPathways_2024_Human",
    "MSigDB_Hallmark_2020",
]
TOP_N = 20


def load_ortholog_map(cache_path: Path) -> dict[str, str]:
    with open(cache_path) as fh:
        raw = json.load(fh)
    mapping: dict[str, str] = {}
    for key, val in raw.items():
        if val is None:
            continue
        rodent_sym = key.split(":")[0]
        human_sym  = val.get("human_symbol", "")
        if human_sym:
            mapping[rodent_sym] = human_sym
    return mapping


def to_human(rodent_symbols: list[str], orth_map: dict[str, str]) -> list[str]:
    seen: set[str] = set()
    out:  list[str] = []
    for sym in rodent_symbols:
        h = orth_map.get(sym)
        if h and h not in seen:
            seen.add(h)
            out.append(h)
    return out


def enrichr_run(gene_list: list[str], library: str,
                list_name: str = "query") -> pd.DataFrame | None:
    if not gene_list:
        return None
    files = {"list": (None, "\n".join(gene_list)), "description": (None, list_name)}
    try:
        resp = requests.post(f"{ENRICHR_URL}/addList", files=files, timeout=60)
        resp.raise_for_status()
        user_list_id = resp.json()["userListId"]
    except Exception as exc:
        logger.error("  addList failed: %s", exc)
        return None

    try:
        resp = requests.get(
            f"{ENRICHR_URL}/enrich?userListId={user_list_id}&backgroundType={library}",
            timeout=120,
        )
        resp.raise_for_status()
        data = resp.json().get(library, [])
    except Exception as exc:
        logger.error("  enrich failed: %s", exc)
        return None

    if not data:
        return pd.DataFrame()

    cols = ["rank","term","pval","z_score","combined_score","overlapping_genes",
            "adjusted_pval","old_pval","old_adjusted_pval"]
    df = pd.DataFrame(data, columns=cols)
    df["library"]   = library
    df["gene_list"] = list_name
    df["n_genes"]   = df["overlapping_genes"].apply(len)
    return df


def dot_plot(df: pd.DataFrame, title: str, out_path: Path, top_n: int = TOP_N) -> None:
    sub = df[df["adjusted_pval"] > 0].nsmallest(top_n, "adjusted_pval").copy()
    if sub.empty:
        return
    sub["-log10p"] = -np.log10(sub["adjusted_pval"])
    sub = sub.sort_values("-log10p")

    fig, ax = plt.subplots(figsize=(8, max(4, len(sub) * 0.32 + 1.5)))
    sc = ax.scatter(
        sub["-log10p"], range(len(sub)),
        s=sub["n_genes"] * 8,
        c=sub["combined_score"],
        cmap="RdYlBu_r", alpha=0.85, edgecolors="grey", linewidths=0.4,
    )
    plt.colorbar(sc, ax=ax, label="Combined score")
    ax.set_yticks(range(len(sub)))
    labels = [t[:55] + "…" if len(t) > 55 else t for t in sub["term"]]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log₁₀(adjusted p-value)", fontsize=10)
    ax.axvline(-np.log10(0.05), ls="--", lw=0.8, color="grey", alpha=0.6)
    ax.set_title(title, fontsize=10, pad=8)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("  Saved %s", out_path.name)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    orth_map = load_ortholog_map(CACHE_JSON)
    logger.info("Ortholog map: %d entries", len(orth_map))

    models = ["CCI", "SNI", "SNL", "CFA"]
    all_results: dict[str, list[pd.DataFrame]] = {}

    for model in models:
        f = STRAT_DIR / f"{model}_pooled.csv"
        if not f.exists():
            logger.warning("Missing: %s", f.name)
            continue
        df = pd.read_csv(f)
        sig = df[(df["padj"] < PADJ_THRESH) & (df["k"] >= K_MIN)]
        genes_rodent = sig.sort_values("pval")["feature_id"].tolist()
        genes_human  = to_human(genes_rodent, orth_map)

        logger.info("\n=== %s: %d sig genes (k>=%d, padj<%.2f) → %d human HGNC ===",
                    model, len(genes_rodent), K_MIN, PADJ_THRESH, len(genes_human))

        if len(genes_human) < MIN_LIST:
            logger.info("  Too few genes for Enrichr ORA — reporting descriptively only")
            print(f"\n[{model}] {len(genes_rodent)} sig genes (below ORA threshold):")
            for g in genes_rodent[:10]:
                row = sig[sig["feature_id"] == g].iloc[0]
                print(f"  {g}: yi={row['yi']:+.3f}, k={row['k']}, padj={row['padj']:.2e}")
            continue

        all_results[model] = []
        for lib in LIBRARIES:
            logger.info("  Enrichr: %s", lib)
            res = enrichr_run(genes_human, lib, list_name=f"{model}_{lib}")
            if res is None or res.empty:
                time.sleep(0.5)
                continue
            sig_count = (res["adjusted_pval"] < 0.05).sum()
            logger.info("    %d terms (padj<0.05: %d)", len(res), sig_count)

            res.to_csv(OUT_DIR / f"{model}_{lib}.csv", index=False)
            dot_plot(
                res,
                title=f"{model} — {lib}  ({len(genes_human)} genes)",
                out_path=OUT_DIR / f"{model}_{lib}_top{TOP_N}.pdf",
            )
            all_results[model].append(res)
            time.sleep(0.5)

    # ---- Model comparison table (top 5 per model × library) ----------------
    comparison_rows = []
    for model, result_list in all_results.items():
        for res in result_list:
            top5 = res[res["adjusted_pval"] < 0.05].nsmallest(5, "adjusted_pval")
            for _, row in top5.iterrows():
                comparison_rows.append({
                    "model":         model,
                    "library":       row["library"],
                    "term":          row["term"],
                    "adjusted_pval": row["adjusted_pval"],
                    "n_genes":       row["n_genes"],
                    "combined_score":row["combined_score"],
                })

    if comparison_rows:
        comp_df = pd.DataFrame(comparison_rows).sort_values(["model", "adjusted_pval"])
        comp_df.to_csv(OUT_DIR / "model_comparison.csv", index=False)
        logger.info("\nModel comparison table → model_comparison.csv (%d rows)", len(comp_df))

    # ---- Summary printout --------------------------------------------------
    print("\n========== Per-Model Pathway Enrichment Summary ==========")
    for model, result_list in all_results.items():
        combined = pd.concat(result_list, ignore_index=True) if result_list else pd.DataFrame()
        if combined.empty:
            print(f"\n[{model}] No significant enrichment")
            continue
        sig = combined[combined["adjusted_pval"] < 0.05]
        print(f"\n[{model}] {len(sig)} significant terms (padj<0.05) "
              f"across {len(LIBRARIES)} libraries")
        top3 = sig.nsmallest(3, "adjusted_pval")[["library","term","adjusted_pval","n_genes"]]
        for _, row in top3.iterrows():
            print(f"  [{row['library'][:8]}] {row['term'][:60]}  "
                  f"padj={row['adjusted_pval']:.2e}  n={row['n_genes']}")

    logger.info("\nStep 4.2 complete. Output: %s", OUT_DIR)


if __name__ == "__main__":
    main()
