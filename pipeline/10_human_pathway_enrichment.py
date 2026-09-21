"""Three-way pathway comparison: rodent-only, human-only, and cross-species concordant.

This is the translational core of the paper (Step 1.8 / Phase 4).

Inputs:
  results/meta/transcriptomics/pooled_effects.csv           — rodent meta-analysis
  results/meta/transcriptomics/stratified/*_human_pooled.csv — human stratified results
  results/cross_species/transcriptomics/concordance_*.csv   — cross-species concordance
  data/interim/ortholog_cache.json                          — rodent → human ortholog map

Strategy:
  1. Rodent list  : k >= K_MIN_ROD, pval < PVAL_ROD (same as step 09 "strong" list)
  2. Human list   : any gene with padj < PADJ_HUMAN in >= 1 human model (gene symbols only;
                    Ensembl IDs from GSE250152 are batch-mapped via mygene.info)
  3. Concordant   : genes that are both_significant=True AND concordant=True in either
                    musculus or norvegicus concordance table
  All three lists use human HGNC symbols for Enrichr.

Outputs:
  results/pathway_enrichment/human/
    gene_list_human.csv           — human-only gene list (HGNC)
    gene_list_concordant.csv      — concordant gene list
    human_{library}.csv
    human_{library}_top20.pdf
    concordant_{library}.csv      (if ≥ 5 genes)
    threeway_comparison.csv       — pathway presence across all three arms
    threeway_venn_top20.pdf       — dot-plot showing overlap across arms
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

REPO_ROOT     = Path(__file__).resolve().parent.parent
POOLED_CSV    = REPO_ROOT / "results" / "meta" / "transcriptomics" / "pooled_effects.csv"
STRAT_DIR     = REPO_ROOT / "results" / "meta" / "transcriptomics" / "stratified"
CROSS_DIR     = REPO_ROOT / "results" / "cross_species" / "transcriptomics"
CACHE_JSON    = REPO_ROOT / "data" / "interim" / "ortholog_cache.json"
OUT_DIR       = REPO_ROOT / "results" / "pathway_enrichment" / "human"

K_MIN_ROD     = 5
PVAL_ROD      = 0.01
PADJ_HUMAN    = 0.05
PVAL_HUMAN    = 0.01    # fallback for underpowered models (lbp_human)

HUMAN_MODELS  = ["neuropathic_human", "lbp_human", "nociplastic_human"]

ENRICHR_URL   = "https://maayanlab.cloud/Enrichr"
MYGENE_URL    = "https://mygene.info/v3/gene"
LIBRARIES     = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "WikiPathways_2024_Human",
    "MSigDB_Hallmark_2020",
]
TOP_N         = 20
BATCH_SIZE    = 1000   # mygene.info batch limit


# ---------------------------------------------------------------------------
# Ensembl → HGNC conversion via mygene.info
# ---------------------------------------------------------------------------

def ensembl_to_hgnc(ensembl_ids: list[str]) -> dict[str, str]:
    """Batch-convert Ensembl gene IDs to HGNC symbols using mygene.info."""
    mapping: dict[str, str] = {}
    for start in range(0, len(ensembl_ids), BATCH_SIZE):
        batch = ensembl_ids[start : start + BATCH_SIZE]
        try:
            resp = requests.post(
                MYGENE_URL,
                data={"ids": ",".join(batch), "fields": "symbol", "species": "human"},
                timeout=30,
            )
            resp.raise_for_status()
            for hit in resp.json():
                q = hit.get("query", "")
                sym = hit.get("symbol", "")
                if q and sym and not hit.get("notfound"):
                    mapping[q] = sym
        except Exception as exc:
            logger.warning("mygene.info batch failed: %s", exc)
        time.sleep(0.3)
    logger.info("  Ensembl→HGNC: %d/%d mapped", len(mapping), len(ensembl_ids))
    return mapping


# ---------------------------------------------------------------------------
# Ortholog map (rodent → human)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Build human gene list from stratified meta-analysis
# ---------------------------------------------------------------------------

def build_human_gene_list() -> list[str]:
    """Return unique HGNC symbols from all human stratified meta-analyses."""
    all_syms: set[str] = set()
    ensembl_ids: list[str] = []

    for model in HUMAN_MODELS:
        f = STRAT_DIR / f"{model}_pooled.csv"
        if not f.exists():
            logger.warning("  Missing: %s", f.name)
            continue
        df = pd.read_csv(f)

        # Choose threshold: use padj if enough genes, else pval fallback
        if "padj" in df.columns:
            sig = df[df["padj"] < PADJ_HUMAN]
        else:
            sig = df[df["pval"] < PVAL_HUMAN]

        # If padj gives <10 genes, fall back to pval
        if len(sig) < 10 and "pval" in df.columns:
            sig = df[df["pval"] < PVAL_HUMAN]
            logger.info("  %s: using pval<%.2f (%d genes)", model, PVAL_HUMAN, len(sig))
        else:
            logger.info("  %s: %d sig genes", model, len(sig))

        for fid in sig["feature_id"]:
            if str(fid).startswith("ENSG"):
                ensembl_ids.append(fid)
            else:
                all_syms.add(str(fid))

    # Map Ensembl IDs to HGNC
    ensembl_ids = list(set(ensembl_ids))
    if ensembl_ids:
        logger.info("  Mapping %d Ensembl IDs via mygene.info …", len(ensembl_ids))
        ens_map = ensembl_to_hgnc(ensembl_ids)
        all_syms.update(ens_map.values())

    result = sorted(all_syms)
    logger.info("Human gene list: %d HGNC symbols total", len(result))
    return result


# ---------------------------------------------------------------------------
# Build concordant gene list
# ---------------------------------------------------------------------------

def build_concordant_list() -> list[str]:
    """Return HGNC symbols concordant and both-significant in mouse or rat."""
    syms: set[str] = set()
    for fname in ["concordance_musculus.csv", "concordance_norvegicus.csv"]:
        f = CROSS_DIR / fname
        if not f.exists():
            continue
        df = pd.read_csv(f)
        mask = df["both_significant"] & df["concordant"]
        syms.update(df.loc[mask, "human_feature_id"].tolist())
    result = sorted(syms)
    logger.info("Concordant gene list: %d genes", len(result))
    return result


# ---------------------------------------------------------------------------
# Enrichr ORA
# ---------------------------------------------------------------------------

def enrichr_run(gene_list: list[str], library: str,
                list_name: str = "query") -> pd.DataFrame | None:
    if not gene_list:
        logger.warning("  Empty gene list for %s / %s — skipping", list_name, library)
        return None

    files = {"list": (None, "\n".join(gene_list)), "description": (None, list_name)}
    try:
        resp = requests.post(f"{ENRICHR_URL}/addList", files=files, timeout=60)
        resp.raise_for_status()
        user_list_id = resp.json()["userListId"]
    except Exception as exc:
        logger.error("  Enrichr addList failed: %s", exc)
        return None

    try:
        resp = requests.get(
            f"{ENRICHR_URL}/enrich?userListId={user_list_id}&backgroundType={library}",
            timeout=120,
        )
        resp.raise_for_status()
        data = resp.json().get(library, [])
    except Exception as exc:
        logger.error("  Enrichr enrich failed: %s", exc)
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


# ---------------------------------------------------------------------------
# Dot plot
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Three-way comparison table
# ---------------------------------------------------------------------------

def build_threeway_comparison(
    rodent_results: list[pd.DataFrame],
    human_results:  list[pd.DataFrame],
    concordant_results: list[pd.DataFrame],
    out_path: Path,
) -> pd.DataFrame:
    """Merge Enrichr results across three arms and flag pathway presence."""
    def sig_terms(dfs: list[pd.DataFrame]) -> pd.DataFrame:
        if not dfs:
            return pd.DataFrame(columns=["library", "term", "adjusted_pval",
                                         "combined_score", "n_genes",
                                         "overlapping_genes"])
        combined = pd.concat(dfs, ignore_index=True)
        return combined[combined["adjusted_pval"] < 0.05]

    rod_sig  = sig_terms(rodent_results)
    hum_sig  = sig_terms(human_results)
    con_sig  = sig_terms(concordant_results)

    rod_key = set(zip(rod_sig["library"], rod_sig["term"])) if not rod_sig.empty else set()
    hum_key = set(zip(hum_sig["library"], hum_sig["term"])) if not hum_sig.empty else set()
    con_key = set(zip(con_sig["library"], con_sig["term"])) if not con_sig.empty else set()

    all_keys = rod_key | hum_key | con_key
    rows = []
    for lib, term in sorted(all_keys):
        def best_padj(df_sig: pd.DataFrame) -> float:
            sub = df_sig[(df_sig["library"] == lib) & (df_sig["term"] == term)]
            return float(sub["adjusted_pval"].min()) if not sub.empty else np.nan

        rows.append({
            "library":          lib,
            "term":             term,
            "in_rodent":        (lib, term) in rod_key,
            "in_human":         (lib, term) in hum_key,
            "in_concordant":    (lib, term) in con_key,
            "padj_rodent":      best_padj(rod_sig),
            "padj_human":       best_padj(hum_sig),
            "padj_concordant":  best_padj(con_sig),
        })

    df_out = pd.DataFrame(rows)
    # Classify into arms
    def arm(row: pd.Series) -> str:
        arms = []
        if row["in_rodent"]:
            arms.append("rodent")
        if row["in_human"]:
            arms.append("human")
        if row["in_concordant"]:
            arms.append("concordant")
        return "+".join(arms) if arms else "none"

    df_out["arm"] = df_out.apply(arm, axis=1)
    df_out = df_out.sort_values(["arm", "padj_human", "padj_rodent"])
    df_out.to_csv(out_path, index=False)
    logger.info("Three-way comparison: %d terms saved → %s", len(df_out), out_path.name)
    return df_out


# ---------------------------------------------------------------------------
# Three-way dot plot
# ---------------------------------------------------------------------------

def threeway_dot_plot(df: pd.DataFrame, out_path: Path, top_n: int = 20) -> None:
    """Dot plot showing top pathways colored by arm membership."""
    if df.empty:
        return

    arm_order = ["rodent+human+concordant", "rodent+human", "human+concordant",
                 "rodent+concordant", "human", "rodent", "concordant"]
    arm_colors = {
        "rodent+human+concordant": "#7b2d8b",
        "rodent+human":            "#2166ac",
        "human+concordant":        "#d6604d",
        "rodent+concordant":       "#4dac26",
        "human":                   "#f4a582",
        "rodent":                  "#92c5de",
        "concordant":              "#b8e186",
    }

    # Pick top_n terms: prioritise multi-arm, then smallest padj
    df = df.copy()
    df["_arm_order"] = df["arm"].map({a: i for i, a in enumerate(arm_order)}).fillna(99)
    df["_best_padj"] = df[["padj_rodent","padj_human","padj_concordant"]].min(axis=1)
    top = df.sort_values(["_arm_order","_best_padj"]).head(top_n).copy()
    top["-log10p"] = -np.log10(top["_best_padj"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(9, max(5, len(top) * 0.35 + 2)))
    colors = [arm_colors.get(a, "grey") for a in top["arm"]]
    ax.scatter(top["-log10p"], range(len(top)), c=colors, s=80, zorder=3)
    ax.set_yticks(range(len(top)))
    labels = [f"[{row['library'][:8]}] {row['term'][:48]}" for _, row in top.iterrows()]
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("-log₁₀(best adjusted p-value)", fontsize=9)
    ax.axvline(-np.log10(0.05), ls="--", lw=0.8, color="grey", alpha=0.6)
    ax.set_title("Three-way pathway comparison\n(rodent / human / concordant)", fontsize=10)

    # Legend
    handles = [
        plt.Line2D([0],[0], marker="o", color="w",
                   markerfacecolor=c, markersize=8, label=a)
        for a, c in arm_colors.items() if a in top["arm"].values
    ]
    ax.legend(handles=handles, fontsize=7, loc="lower right", title="Arm")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Three-way dot plot → %s", out_path.name)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    orth_map = load_ortholog_map(CACHE_JSON)
    logger.info("Ortholog map: %d entries", len(orth_map))

    # ---- Rodent gene list (same criteria as step 09 "strong") ---------------
    logger.info("\n=== Rodent gene list ===")
    rod_df = pd.read_csv(POOLED_CSV)
    sig_rod = rod_df[rod_df["pval_pooled"] < PVAL_ROD]
    strong_rod_rodent = (sig_rod[sig_rod["k"] >= K_MIN_ROD]
                         .sort_values("pval_pooled")["feature_id"].tolist())
    strong_rod_human  = to_human(strong_rod_rodent, orth_map)
    logger.info("Rodent: %d features (k>=%d, pval<%.2f) → %d human HGNC",
                len(strong_rod_rodent), K_MIN_ROD, PVAL_ROD, len(strong_rod_human))

    # ---- Human gene list ----------------------------------------------------
    logger.info("\n=== Human gene list ===")
    human_genes = build_human_gene_list()
    pd.DataFrame({"hgnc": human_genes}).to_csv(OUT_DIR / "gene_list_human.csv", index=False)

    # ---- Concordant gene list -----------------------------------------------
    logger.info("\n=== Concordant gene list ===")
    concordant_genes = build_concordant_list()
    pd.DataFrame({"hgnc": concordant_genes}).to_csv(
        OUT_DIR / "gene_list_concordant.csv", index=False)

    # ---- Enrichr ORA --------------------------------------------------------
    rodent_results:     list[pd.DataFrame] = []
    human_results:      list[pd.DataFrame] = []
    concordant_results: list[pd.DataFrame] = []

    for list_name, gene_list, result_bucket, out_prefix in [
        ("rodent_strong", strong_rod_human,  rodent_results,     "rodent"),
        ("human_all",     human_genes,        human_results,      "human"),
        ("concordant",    concordant_genes,   concordant_results, "concordant"),
    ]:
        if len(gene_list) < 5:
            logger.info("\n=== %s: only %d genes — skipping Enrichr ===", list_name, len(gene_list))
            continue
        logger.info("\n===== Enrichr: %s (%d genes) =====", list_name, len(gene_list))
        for lib in LIBRARIES:
            logger.info("  %s", lib)
            res = enrichr_run(gene_list, lib, list_name=f"{list_name}_{lib}")
            if res is None or res.empty:
                time.sleep(0.5)
                continue
            sig_count = (res["adjusted_pval"] < 0.05).sum()
            logger.info("    %d terms (padj<0.05: %d)", len(res), sig_count)

            out_csv = OUT_DIR / f"{out_prefix}_{lib}.csv"
            res.to_csv(out_csv, index=False)

            dot_plot(
                res,
                title=f"{list_name} — {lib}",
                out_path=OUT_DIR / f"{out_prefix}_{lib}_top{TOP_N}.pdf",
            )
            result_bucket.append(res)
            time.sleep(0.5)

    # ---- Three-way comparison -----------------------------------------------
    logger.info("\n=== Three-way comparison ===")
    threeway = build_threeway_comparison(
        rodent_results, human_results, concordant_results,
        out_path=OUT_DIR / "threeway_comparison.csv",
    )

    if not threeway.empty:
        threeway_dot_plot(
            threeway,
            out_path=OUT_DIR / "threeway_venn_top20.pdf",
        )

    # ---- Summary printout ---------------------------------------------------
    print("\n========== Human Pathway Enrichment Summary ==========")
    print("\nGene lists:")
    print(f"  Rodent (strong, k>={K_MIN_ROD}, pval<{PVAL_ROD}): {len(strong_rod_human)} HGNC")
    print(f"  Human (all models, padj<{PADJ_HUMAN} or pval<{PVAL_HUMAN}): {len(human_genes)} HGNC")
    print(f"  Concordant (both-sig + concordant): {len(concordant_genes)} genes")

    if not threeway.empty:
        arm_counts = threeway.groupby("arm").size()
        print(f"\nThree-way pathway overlap ({len(threeway)} total sig terms):")
        for arm, n in arm_counts.sort_index().items():
            print(f"  {arm:30s} {n} pathways")

        shared = threeway[threeway["in_rodent"] & threeway["in_human"]]
        if not shared.empty:
            print(f"\nShared rodent+human pathways (n={len(shared)}):")
            for _, row in shared.head(10).iterrows():
                print(f"  [{row['library'][:10]}] {row['term'][:55]}")
                print(f"    padj_rodent={row['padj_rodent']:.2e}  "
                      f"padj_human={row['padj_human']:.2e}")

    logger.info("\nStep 10 complete. Output: %s", OUT_DIR)


if __name__ == "__main__":
    main()
