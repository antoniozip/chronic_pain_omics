"""Cross-species integration pipeline step.

Reads pooled meta-analysis results (results/meta/{modality}/pooled_effects.csv),
splits into human-only and animal-only study pools using the harmonized manifest,
maps animal gene symbols to human orthologs via MyGene.info (cached), and
computes concordance metrics.

Output:
    results/cross_species/{modality}/concordance_{species}.csv
    results/cross_species/{modality}/summary.csv   (correlation + concordance stats)

Usage:
    uv run python pipeline/07_cross_species.py --modality transcriptomics
    uv run python pipeline/07_cross_species.py --modality all
    uv run python pipeline/07_cross_species.py --modality transcriptomics --dry-run
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict
from pathlib import Path

import pandas as pd
import yaml

from cp_multiomics.ortholog import OrthologMapper
from cp_multiomics.ortholog.concordance import build_concordance_table, compute_summary

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
ALL_MODALITIES = ["transcriptomics", "proteomics", "metabolomics", "lipidomics"]
# Genomics (GWAS Catalog) is human-only by definition — excluded from cross-species
ANIMAL_SPECIES = ["Mus musculus", "Rattus norvegicus"]
HUMAN_SPECIES = "Homo sapiens"


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_pooled(meta_dir: Path, modality: str) -> pd.DataFrame | None:
    path = meta_dir / modality / "pooled_effects.csv"
    if not path.exists():
        logger.warning("[%s] No pooled effects at %s — skipping", modality, path)
        return None
    df = pd.read_csv(path)
    logger.info("[%s] Loaded %d pooled features", modality, len(df))
    return df


def ortholog_mapping_stats(
    animal_df: pd.DataFrame, species: str, mapper: OrthologMapper
) -> dict:
    """Coverage of the ortholog map over one species' features.

    Reported alongside the concordance so the figure quoted in the paper is
    always the one the pipeline computed. Coverage is low by construction and
    should not be read as a quality problem: most identifiers in the rodent
    tables are array probes, LOC entries, Riken clones and pseudogenes with no
    human counterpart, not named genes.

    The mapper is cache-backed, so this re-lookup costs no network traffic once
    build_concordance_table has run.
    """
    if animal_df.empty or "feature_id" not in animal_df.columns:
        return {"n_animal_features": 0, "n_mapped_to_human": 0,
                "pct_mapped_to_human": float("nan")}

    symbols = animal_df["feature_id"].dropna().astype(str).unique().tolist()
    records = mapper.map_to_human(symbols, from_species=species)
    mapped = sum(
        1 for s in symbols
        if records.get(s) is not None and getattr(records[s], "human_symbol", None)
    )
    return {
        "n_animal_features": len(symbols),
        "n_mapped_to_human": mapped,
        "pct_mapped_to_human": 100 * mapped / len(symbols) if symbols else float("nan"),
    }


def load_species_pooled(
    meta_dir: Path, modality: str, species: str
) -> pd.DataFrame | None:
    """Load one species' pooled effects from 06g_species_meta.R, if present.

    These are the only species-resolved estimates in the pipeline. The
    whole-cohort table pools every study of a feature regardless of species,
    so mouse and rat studies of the same rodent symbol collapse into one row
    and the two species cannot be told apart from it.

    06g writes the generic meta column names (yi, padj); rename them to the
    pooled schema the concordance builder consumes.
    """
    path = meta_dir / modality / "by_species" / f"{species.replace(' ', '_')}_pooled.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path).rename(
        columns={"yi": "yi_pooled", "se": "se_pooled",
                 "pval": "pval_pooled", "padj": "padj_pooled"}
    )
    logger.info("[%s] %s: %d species-resolved features", modality, species, len(df))
    return df


def split_by_species(
    pooled_df: pd.DataFrame,
    harmonized_dir: Path,
    modality: str,
    target_species: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split pooled effects into human-only and target-species-only subsets.

    Uses study_id column in pooled effects joined with harmonized manifest
    to determine which studies contributed human vs. animal effects. The two
    subsets are disjoint by construction: a feature whose pooled estimate
    draws on both human and animal studies is attributed to neither.
    """
    manifest_path = harmonized_dir / modality / "harmonized.jsonl"
    if not manifest_path.exists():
        logger.warning("No harmonized manifest for %s — assuming all human.", modality)
        return pooled_df, pd.DataFrame()

    with open(manifest_path) as f:
        records = [json.loads(line) for line in f if line.strip()]
    manifest = pd.DataFrame(records)

    human_accessions = set(
        manifest.loc[manifest["has_human"] & ~manifest["has_animal"], "accession"]
    )
    animal_accessions = set(
        manifest.loc[
            manifest["has_animal"] &
            manifest["species_canonical"].apply(
                lambda s: target_species in (s if isinstance(s, list) else [])
            ),
            "accession"
        ]
    )
    # Every animal study, not just the target species: a row pooling human and
    # rat studies is not a human estimate merely because the target is mouse.
    any_animal_accessions = set(manifest.loc[manifest["has_animal"], "accession"])

    # pooled_effects has a 'k' column but not per-study ids; we can only split
    # if the meta-analysis output was stratified. Fallback: use all features.
    if "study_ids" in pooled_df.columns:
        # A pooled row is one estimate over all its contributing studies, with
        # no species breakdown. Testing only for intersection puts a row that
        # pools human *and* animal studies on both sides, so the concordance
        # join compares that estimate with itself — concordant by construction,
        # and it is the both-significant pairs (the only ones ever named) that
        # such rows contaminate. Attribute a row to a species only when every
        # contributing study is of that species; a genuinely cross-species row
        # is not evidence about either and is dropped from both subsets.
        study_sets = pooled_df["study_ids"].apply(lambda ids: set(str(ids).split(";")))
        human_mask = study_sets.apply(
            lambda s: bool(s & human_accessions) and not (s & any_animal_accessions)
        )
        animal_mask = study_sets.apply(
            lambda s: bool(s & animal_accessions) and not (s & human_accessions)
        )
        mixed = study_sets.apply(
            lambda s: bool(s & human_accessions) and bool(s & any_animal_accessions)
        )
        logger.info(
            "[%s/%s] %d human-only, %d animal-only rows "
            "(%d rows pool both species and are excluded from both sides)",
            modality, target_species, int(human_mask.sum()), int(animal_mask.sum()),
            int(mixed.sum()),
        )
        return pooled_df[human_mask], pooled_df[animal_mask]

    # No per-row study_ids: return full df for human; empty for animal
    # (cross-species comparison requires stratified meta-analysis output)
    logger.warning(
        "[%s] pooled_effects lacks 'study_ids' column — "
        "cross-species split not possible without stratified meta-analysis.",
        modality,
    )
    return pooled_df, pd.DataFrame()


def run_modality(
    modality: str,
    cfg: dict,
    meta_dir: Path,
    harmonized_dir: Path,
    cs_dir: Path,
    cache_path: Path,
    dry_run: bool,
) -> None:
    pooled = load_pooled(meta_dir, modality)
    human_species_df = load_species_pooled(meta_dir, modality, HUMAN_SPECIES)

    # The whole-cohort table is only the fallback: with species-resolved pools
    # from 06g this step never reads it. Returning early when it is absent
    # therefore skipped runs that had everything they needed -- which is what
    # every cross-species sensitivity arm did. Each arm re-runs 06g into a work
    # directory holding only by_species/, so 07 found no pooled_effects.csv
    # there, logged one warning, exited 0, and the driver's copy step then
    # found nothing to copy and carried on. Table 4's two sensitivity arms had
    # been stale since 2026-08-29 for that reason, through two chains that
    # reported success.
    if pooled is None or pooled.empty:
        if human_species_df is None:
            logger.warning(
                "[%s] neither pooled_effects.csv nor a human by_species pool "
                "— nothing to compare against; skipping", modality)
            return
        logger.info(
            "[%s] no whole-cohort pooled table; using the species-resolved "
            "pools alone, which is what the comparison uses anyway", modality)

    if dry_run:
        logger.info("[dry-run][%s] Would map orthologs for %d features", modality,
                    0 if pooled is None else len(pooled))
        return

    mapper = OrthologMapper(cache_path)
    summaries: list[dict] = []

    for species in ANIMAL_SPECIES:
        logger.info("[%s] Processing species: %s", modality, species)

        animal_species_df = load_species_pooled(meta_dir, modality, species)
        if human_species_df is not None and animal_species_df is not None:
            human_df, animal_df = human_species_df, animal_species_df
        else:
            # Falls back to carving the whole-cohort table by contributing
            # study. That cannot separate mouse from rat, since both pool into
            # the same rodent-symbol row, so the per-species results are only
            # as distinct as the study sets behind them.
            logger.warning(
                "[%s/%s] No species-stratified pooled effects (run "
                "pipeline/06g_species_meta.R) — falling back to splitting the "
                "whole-cohort table, which cannot separate rodent species.",
                modality, species,
            )
            if pooled is None:
                logger.warning(
                    "[%s/%s] no species-resolved pool and no whole-cohort "
                    "table to carve; skipping", modality, species)
                continue
            human_df, animal_df = split_by_species(
                pooled, harmonized_dir, modality, species
            )

        if animal_df.empty:
            logger.info("[%s] No animal-only studies for %s — skipping.", modality, species)
            continue

        concordance = build_concordance_table(
            human_df=human_df,
            animal_df=animal_df,
            animal_species=species,
            mapper=mapper,
            padj_threshold=cfg.get("da", {}).get("padj_threshold", 0.05),
        )

        if concordance.empty:
            continue

        out_dir = cs_dir / modality
        out_dir.mkdir(parents=True, exist_ok=True)
        species_slug = species.split()[1].lower()
        concordance.to_csv(out_dir / f"concordance_{species_slug}.csv", index=False)
        logger.info("[%s/%s] %d ortholog pairs written.", modality, species, len(concordance))

        mapping = ortholog_mapping_stats(animal_df, species, mapper)
        logger.info(
            "[%s/%s] ortholog coverage: %d/%d features map to human (%.1f%%)",
            modality, species, mapping["n_mapped_to_human"],
            mapping["n_animal_features"], mapping["pct_mapped_to_human"],
        )

        summary = compute_summary(concordance, modality)
        summaries.append({**asdict(summary), "species": species, **mapping})

    if summaries:
        pd.DataFrame(summaries).to_csv(cs_dir / modality / "summary.csv", index=False)
        logger.info("[%s] Summary written.", modality)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--modality", choices=ALL_MODALITIES + ["all"], required=True)
    p.add_argument("--analysis-config", type=Path,
                   default=REPO_ROOT / "conf" / "analysis" / "default.yaml")
    p.add_argument("--meta-dir",       type=Path, default=REPO_ROOT / "results" / "meta")
    p.add_argument("--harmonized-dir", type=Path, default=REPO_ROOT / "data" / "interim")
    p.add_argument("--out-dir",        type=Path, default=REPO_ROOT / "results" / "cross_species")
    p.add_argument("--cache-path",     type=Path,
                   default=REPO_ROOT / "data" / "interim" / "ortholog_cache.json")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg  = load_config(args.analysis_config)
    mods = ALL_MODALITIES if args.modality == "all" else [args.modality]

    for mod in mods:
        run_modality(
            modality=mod, cfg=cfg,
            meta_dir=args.meta_dir,
            harmonized_dir=args.harmonized_dir,
            cs_dir=args.out_dir,
            cache_path=args.cache_path,
            dry_run=args.dry_run,
        )

    logger.info("Step 07 complete.")


if __name__ == "__main__":
    main()
