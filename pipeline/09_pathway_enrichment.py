"""Pathway enrichment analysis on multi-study transcriptomics features.

Inputs:
  results/meta/transcriptomics/pooled_effects.csv
  data/interim/ortholog_cache.json

Strategy:
  1. Convert rodent gene symbols to human HGNC via the ortholog cache (step 07).
  2. Define three ranked/filtered gene lists:
       - "strong"    : k >= K_MIN and padj < PADJ_THRESH (primary)
       - "moderate"  : k >= 2 and padj < PADJ_THRESH (extended)
       - "ranked_bg" : all k >= 2 features, rank-ordered by -log10(padj)*sign
  3. Run Enrichr ORA (addList → enrich) for "strong" and "moderate" lists against
     five pathway libraries (KEGG_2021_Human, GO_Biological_Process_2023,
     Reactome_2022, WikiPathways_2024_Human, MSigDB_Hallmark_2020).
  4. Save enrichment tables to results/pathway_enrichment/.
  5. Generate dot-plot figures (top 20 terms per library per gene list).

Outputs:
  results/pathway_enrichment/
    {list}_{library}.csv          — full enrichment results
    {list}_{library}_top20.pdf    — dot plot
    summary_top10.csv             — top 10 hits across libraries, both lists
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

PER_STUDY_CSV  = (Path(__file__).resolve().parent.parent / "results" / "per_study"
                  / "transcriptomics" / "effects_normalized_all.csv")

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

REPO_ROOT      = Path(__file__).resolve().parent.parent
POOLED_CSV     = REPO_ROOT / "results" / "meta" / "transcriptomics" / "pooled_effects.csv"
CACHE_JSON     = REPO_ROOT / "data" / "interim" / "ortholog_cache.json"
OUT_DIR        = REPO_ROOT / "results" / "pathway_enrichment"

K_MIN          = 5      # minimum study count for "strong" list
K_MOD          = 2      # minimum study count for "moderate" list
PVAL_THRESH    = 0.01   # nominal p-value threshold (padj is too stringent given ~47k tests)

ENRICHR_URL    = "https://maayanlab.cloud/Enrichr"
LIBRARIES      = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "WikiPathways_2024_Human",
    "MSigDB_Hallmark_2020",
]
TOP_N          = 20     # terms shown in dot plot


# ---------------------------------------------------------------------------
# Ortholog conversion (rodent → human)
# ---------------------------------------------------------------------------

def load_ortholog_map(cache_path: Path) -> dict[str, str]:
    """Return {rodent_symbol: human_HGNC} from the ortholog cache."""
    with open(cache_path) as fh:
        raw = json.load(fh)
    mapping: dict[str, str] = {}
    for key, val in raw.items():
        if val is None:
            continue
        # cache key is "Symbol:taxid"
        rodent_sym = key.split(":")[0]
        human_sym  = val.get("human_symbol", "")
        if human_sym:
            mapping[rodent_sym] = human_sym
    return mapping


def to_human(rodent_symbols: list[str], orth_map: dict[str, str]) -> list[str]:
    """Map list of rodent symbols to unique human HGNC symbols (order preserved)."""
    seen: set[str] = set()
    out: list[str] = []
    for sym in rodent_symbols:
        h = orth_map.get(sym)
        if h and h not in seen:
            seen.add(h)
            out.append(h)
    return out


# ---------------------------------------------------------------------------
# Enrichr API
# ---------------------------------------------------------------------------

def enrichr_run(gene_list: list[str], library: str,
                list_name: str = "query") -> pd.DataFrame | None:
    """Submit gene list to Enrichr and retrieve enrichment results.

    Returns a DataFrame or None on failure.
    """
    if not gene_list:
        logger.warning("Empty gene list for %s / %s — skipping", list_name, library)
        return None

    # Step 1: upload gene list (Enrichr requires multipart/form-data)
    add_url = f"{ENRICHR_URL}/addList"
    files    = {"list": (None, "\n".join(gene_list)), "description": (None, list_name)}
    try:
        resp = requests.post(add_url, files=files, timeout=60)
        resp.raise_for_status()
        user_list_id = resp.json()["userListId"]
    except Exception as exc:
        logger.error("Enrichr addList failed for %s / %s: %s", list_name, library, exc)
        return None

    # Step 2: retrieve results
    enrich_url = (
        f"{ENRICHR_URL}/enrich?userListId={user_list_id}"
        f"&backgroundType={library}"
    )
    try:
        resp = requests.get(enrich_url, timeout=120)
        resp.raise_for_status()
        data = resp.json().get(library, [])
    except Exception as exc:
        logger.error("Enrichr enrich failed for %s / %s: %s", list_name, library, exc)
        return None

    if not data:
        return pd.DataFrame()

    columns = [
        "rank", "term", "pval", "z_score", "combined_score",
        "overlapping_genes", "adjusted_pval", "old_pval", "old_adjusted_pval",
    ]
    df = pd.DataFrame(data, columns=columns)
    df["library"]   = library
    df["gene_list"] = list_name
    df["n_genes"]   = df["overlapping_genes"].apply(len)
    return df


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def dot_plot(df: pd.DataFrame, title: str, out_path: Path, top_n: int = TOP_N) -> None:
    """Dot plot: x = -log10(adjusted_pval), dot size = n_genes."""
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
        cmap="RdYlBu_r",
        alpha=0.85,
        edgecolors="grey",
        linewidths=0.4,
    )
    plt.colorbar(sc, ax=ax, label="Combined score")
    ax.set_yticks(range(len(sub)))
    # Wrap long term names
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
# Pathway forest plots
# ---------------------------------------------------------------------------

def pathway_forest_plots(
    combined: pd.DataFrame,
    pooled_df: pd.DataFrame,
    orth_map: dict[str, str],
    out_dir: Path,
    top_pathways: int = 6,
) -> None:
    """For the top significant pathways, plot per-study effects of the single
    best-ranked gene contributing to that pathway."""
    if not PER_STUDY_CSV.exists():
        logger.warning("Per-study effects CSV not found — skipping pathway forest plots")
        return

    per_study = pd.read_csv(PER_STUDY_CSV)
    # reverse orth_map to get human→rodent(s)
    human_to_rodent: dict[str, list[str]] = {}
    for rod, hum in orth_map.items():
        human_to_rodent.setdefault(hum, []).append(rod)

    sig = combined[combined["adjusted_pval"] < 0.05].copy()
    if sig.empty:
        return

    # Pick one representative term per library × gene_list combo, avoiding duplicates
    top_terms = (
        sig.sort_values("adjusted_pval")
        .drop_duplicates(subset=["library"])
        .head(top_pathways)
    )

    pdf_path = out_dir / "pathway_representative_forest.pdf"
    import matplotlib.backends.backend_pdf as mpdf
    with mpdf.PdfPages(pdf_path) as pdf:
        for _, term_row in top_terms.iterrows():
            term_name  = term_row["term"]
            lib        = term_row["library"]
            genes_human = [g.strip() for g in term_row["overlapping_genes"] if g.strip()]
            if not genes_human:
                continue

            # Map human → rodent, find the best-replicated gene in per_study data
            rodent_candidates: list[str] = []
            for hg in genes_human:
                rodent_candidates.extend(human_to_rodent.get(hg, [hg]))

            counts = per_study[per_study["feature_id"].isin(rodent_candidates)]\
                .groupby("feature_id")["study_id"].nunique()
            if counts.empty:
                continue
            best_gene = counts.idxmax()
            gene_data = per_study[per_study["feature_id"] == best_gene].copy()
            if len(gene_data) < 2:
                continue

            # Weighted mean for the pooled estimate line
            w = 1.0 / (gene_data["se"] ** 2)
            pooled_yi = (gene_data["effect_size"] * w).sum() / w.sum()
            pooled_se = np.sqrt(1.0 / w.sum())

            gene_data = gene_data.sort_values("effect_size")
            fig, ax = plt.subplots(figsize=(7, max(4, len(gene_data) * 0.35 + 2)))

            y_pos = np.arange(len(gene_data))
            ax.errorbar(
                gene_data["effect_size"], y_pos,
                xerr=1.96 * gene_data["se"],
                fmt="o", color="#2166ac", ecolor="#92c5de",
                markersize=5, linewidth=1, capsize=3,
            )
            ax.axvline(0, color="black", linewidth=0.7, linestyle="--")
            ax.axvline(pooled_yi, color="#d6604d", linewidth=1.5, linestyle="-",
                       label=f"Pooled: {pooled_yi:+.3f} ± {pooled_se:.3f}")
            ax.axvspan(pooled_yi - 1.96 * pooled_se, pooled_yi + 1.96 * pooled_se,
                       alpha=0.15, color="#d6604d")
            ax.set_yticks(y_pos)
            ax.set_yticklabels(gene_data["study_id"], fontsize=8)
            ax.set_xlabel("log₂ fold change (case vs. control)", fontsize=9)
            ax.set_title(
                f"{best_gene}  •  {term_name[:50]}\n({lib}  padj={term_row['adjusted_pval']:.2e})",
                fontsize=9,
            )
            ax.legend(fontsize=8)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            logger.info("  Pathway forest: %s (%s)", best_gene, term_name[:40])

    logger.info("Pathway representative forest plots → %s", pdf_path.name)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- Load pooled effects ------------------------------------------------
    df = pd.read_csv(POOLED_CSV)
    logger.info("Loaded %d features from pooled effects", len(df))

    # ---- Load ortholog map --------------------------------------------------
    orth_map = load_ortholog_map(CACHE_JSON)
    logger.info("Ortholog map: %d rodent → human entries", len(orth_map))

    # ---- Build gene lists ---------------------------------------------------
    # Use nominal p-value: multiple-testing correction over ~47k features leaves
    # only ~12 genes at padj<0.05, which is too small for ORA. Nominal pval<0.01
    # with replication filters gives interpretable gene lists.
    sig = df[df["pval_pooled"] < PVAL_THRESH].copy()
    strong_rodent   = sig[sig["k"] >= K_MIN].sort_values("pval_pooled")["feature_id"].tolist()
    moderate_rodent = sig[sig["k"] >= K_MOD].sort_values("pval_pooled")["feature_id"].tolist()

    strong_human   = to_human(strong_rodent,   orth_map)
    moderate_human = to_human(moderate_rodent, orth_map)

    logger.info(
        "Gene lists  strong: %d rodent (k>=%d) → %d human  |  moderate: %d (k>=%d) → %d human",
        len(strong_rodent), K_MIN, len(strong_human),
        len(moderate_rodent), K_MOD, len(moderate_human),
    )

    # Save gene lists for reference
    pd.DataFrame({"rodent": strong_rodent}).to_csv(
        OUT_DIR / "gene_list_strong_rodent.csv", index=False)
    pd.DataFrame({"human": strong_human}).to_csv(
        OUT_DIR / "gene_list_strong_human.csv", index=False)
    pd.DataFrame({"rodent": moderate_rodent}).to_csv(
        OUT_DIR / "gene_list_moderate_rodent.csv", index=False)
    pd.DataFrame({"human": moderate_human}).to_csv(
        OUT_DIR / "gene_list_moderate_human.csv", index=False)

    # ---- Enrichr ORA --------------------------------------------------------
    all_results: list[pd.DataFrame] = []

    for gene_list_name, gene_list in [("strong", strong_human), ("moderate", moderate_human)]:
        logger.info("\n===== Gene list: %s (%d genes) =====", gene_list_name, len(gene_list))
        for lib in LIBRARIES:
            logger.info("  Enrichr: %s", lib)
            res = enrichr_run(gene_list, lib, list_name=f"{gene_list_name}_{lib}")
            if res is None or res.empty:
                logger.info("    No results.")
                time.sleep(0.5)
                continue

            sig_res = res[res["adjusted_pval"] < 0.05]
            logger.info("    %d terms (padj<0.05: %d)", len(res), len(sig_res))

            out_csv = OUT_DIR / f"{gene_list_name}_{lib}.csv"
            res.to_csv(out_csv, index=False)

            # Dot plot for top N
            out_pdf = OUT_DIR / f"{gene_list_name}_{lib}_top{TOP_N}.pdf"
            dot_plot(
                res, title=f"{gene_list_name.capitalize()} genes — {lib}",
                out_path=out_pdf, top_n=TOP_N,
            )

            all_results.append(res)
            time.sleep(0.5)   # be polite to the API

    # ---- Summary table -------------------------------------------------------
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        top10 = (
            combined[combined["adjusted_pval"] < 0.05]
            .sort_values(["gene_list", "adjusted_pval"])
            .groupby(["gene_list", "library"])
            .head(10)
            [["gene_list", "library", "rank", "term", "pval",
              "adjusted_pval", "n_genes", "combined_score", "overlapping_genes"]]
        )
        top10.to_csv(OUT_DIR / "summary_top10.csv", index=False)
        logger.info("\nSummary top-10 table saved: %d rows", len(top10))

        # Print quick overview to stdout
        print("\n========== Enrichment Overview ==========")
        for (gl, lib), grp in (
            combined[combined["adjusted_pval"] < 0.05]
            .groupby(["gene_list", "library"])
        ):
            best = grp.nsmallest(3, "adjusted_pval")[["term", "adjusted_pval", "n_genes"]]
            print(f"\n[{gl}] {lib}  ({len(grp)} significant terms)")
            for _, row in best.iterrows():
                print(f"  {row['term'][:60]:60s}  padj={row['adjusted_pval']:.2e}  "
                      f"n={row['n_genes']}")

        # ---- Pathway representative forest plots ----------------------------
        pathway_forest_plots(combined, df, orth_map, OUT_DIR, top_pathways=6)

    logger.info("\nStep 09 complete. Output: %s", OUT_DIR)


if __name__ == "__main__":
    main()
