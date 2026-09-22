"""Paper figure generation pipeline step.

Reads meta-analysis results and cross-species concordance tables to produce
all manuscript figures in manuscript/figures/.

Figures produced:
    {modality}_volcano.pdf          — per-modality volcano plots (step 06 output)
    {modality}_concordance_{sp}.pdf — human vs animal concordance scatter plots
    prisma_flow.pdf                 — PRISMA 2020 flow diagram
    graphical_abstract.pdf          — front-page summary of the three findings

Usage:
    uv run python pipeline/08_figures.py
    uv run python pipeline/08_figures.py --modality transcriptomics
    uv run python pipeline/08_figures.py --prisma-only
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

from cp_multiomics.prisma_flow import metabolomic_units, repository_flows
from cp_multiomics.viz import (
    plot_concordance_scatter,
    plot_graphical_abstract,
    plot_repository_flow,
    plot_volcano,
    set_publication_style,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT    = Path(__file__).parent.parent
ALL_MODALITIES = ["genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics"]
ANIMAL_SPECIES = ["Mus musculus", "Rattus norvegicus"]


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_pooled(meta_dir: Path, modality: str) -> pd.DataFrame:
    path = meta_dir / modality / "pooled_effects.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def load_concordance(cs_dir: Path, modality: str, species: str) -> pd.DataFrame:
    slug = species.split()[1].lower()
    path = cs_dir / modality / f"concordance_{slug}.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def repository_flow_summary(flows: list, metabolomic_units: int | None) -> list[str]:
    """Lines for the closing box: what each arm contributes to the synthesis."""
    by = {f.source: f for f in flows}
    geo, sc = by["GEO (bulk transcriptomics)"], by["GEO single-cell series"]
    gwas, pride = by["GWAS Catalog"], by["PRIDE"]
    metab = by["MetaboLights"].included + by["Metabolomics Workbench"].included
    return [
        f"Transcriptomic: {geo.units} study units from {geo.included} GEO accessions, "
        f"and {sc.included} single-cell series pooled as a separate pseudobulk stratum",
        f"Genomic: {gwas.included} genome-wide association studies",
        f"Proteomic: {pride.units} study units from {pride.included} PRIDE datasets",
        f"Metabolomic: {metabolomic_units} study units from {metab} studies "
        "in both repositories",
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--modality", choices=ALL_MODALITIES + ["all"], default="all")
    p.add_argument("--analysis-config", type=Path,
                   default=REPO_ROOT / "conf" / "analysis" / "default.yaml")
    p.add_argument("--meta-dir",    type=Path, default=REPO_ROOT / "results" / "meta")
    p.add_argument("--cs-dir",      type=Path, default=REPO_ROOT / "results" / "cross_species")
    p.add_argument("--figures-dir", type=Path, default=REPO_ROOT / "manuscript" / "figures")
    p.add_argument("--prisma-only", action="store_true",
                   help="Generate only the PRISMA flow diagram and exit")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg  = load_config(args.analysis_config)
    set_publication_style()

    padj  = cfg.get("da", {}).get("padj_threshold", 0.05)
    lfc   = cfg.get("da", {}).get("lfc_threshold", 0.5)

    # --- PRISMA flow ---
    # Datasets through the repositories that supplied them; see
    # cp_multiomics.prisma_flow for why this is no longer the PubMed search.
    flows = repository_flows(REPO_ROOT)
    plot_repository_flow(
        flows, repository_flow_summary(flows, metabolomic_units(REPO_ROOT)),
        out_dir=args.figures_dir)

    if args.prisma_only:
        return

    # --- Graphical abstract ---
    # Study counts are the corpus sizes the abstract reports; every other
    # number in the figure is read from results/ when it is drawn, so the
    # figure cannot disagree with the meta-analysis it summarises.
    plot_graphical_abstract(
        REPO_ROOT / "results",
        args.figures_dir / "graphical_abstract.pdf",
        arms=[("Transcriptomic", "transcriptomics", 58),
              ("Genomic", "genomics", 57),
              ("Proteomic", "proteomics", 43),
              ("Metabolomic", "metabolomics", 9)],
    )

    modalities = ALL_MODALITIES if args.modality == "all" else [args.modality]

    # --- Per-modality volcano plots ---
    for mod in modalities:
        df = load_pooled(args.meta_dir, mod)
        if df.empty:
            logger.info("[%s] No pooled effects — skipping volcano.", mod)
            continue
        plot_volcano(df, modality=mod, out_dir=args.figures_dir,
                     padj_threshold=padj, lfc_threshold=lfc, top_n_labels=15)

    # --- Cross-species concordance scatter ---
    for mod in [m for m in modalities if m != "genomics"]:
        for species in ANIMAL_SPECIES:
            conc = load_concordance(args.cs_dir, mod, species)
            if conc.empty:
                continue
            plot_concordance_scatter(
                conc, modality=mod, animal_species=species,
                out_dir=args.figures_dir, padj_threshold=padj,
            )

    logger.info("Step 08 complete. Figures written to %s", args.figures_dir)


if __name__ == "__main__":
    main()
