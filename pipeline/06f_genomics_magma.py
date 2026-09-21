"""Run MAGMA gene-based analysis on the pooled chronic-pain METAL meta-analyses.

For each poolable pain group, aggregates the group's SNP-level meta-analysis
(METAL .tbl) to gene-level association with MAGMA, using the 1000G EUR reference
for LD. Gene-based tests recover signal from many sub-threshold SNPs that no single
marker would surface. Emits per-gene results, Bonferroni-significant genes, and a
provenance manifest.

    uv run python pipeline/06f_genomics_magma.py

Complements the SNP-level METAL stage; the manuscript's genomics section gains a
gene-level layer.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

from cp_multiomics.genomics import magma_runner as magma_runner_mod
from cp_multiomics.genomics.magma_runner import (
    annotate,
    parse_genes_out,
    prep_pval_input,
    run_gene_analysis,
    write_snploc,
)
from cp_multiomics.genomics.provenance import write_run_provenance

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "conf" / "genomics" / "magma.yaml"

_MANIFEST_COLUMNS = ["group", "n_snps_in", "n_genes_tested", "n_genes_sig", "status"]


def _default_annotator(config: dict, out_dir: Path):
    """Build the annotation once (cached); return a zero-arg callable -> annot path."""
    def annotator() -> Path:
        annot = out_dir / "g1000_eur.genes.annot"
        if annot.exists():
            logger.info("Using cached annotation %s", annot)
            return annot
        snploc = out_dir / "g1000_eur.snploc"
        write_snploc(Path(config["bfile"] + ".bim"), snploc)
        return annotate(Path(config["magma_bin"]), snploc,
                        Path(config["gene_loc"]), out_dir / "g1000_eur")
    return annotator


def _default_gene_runner(config: dict, annot: Path):
    def gene_runner(group: str, pval_file: Path, out_dir: Path) -> pd.DataFrame | None:
        out_prefix = out_dir / group
        genes_out = run_gene_analysis(
            Path(config["magma_bin"]), Path(config["bfile"]), pval_file,
            annot, out_prefix, n_col="N")
        if not genes_out.exists():
            logger.error("[%s] MAGMA produced no .genes.out", group)
            return None
        return parse_genes_out(genes_out, group)
    return gene_runner


def run(config: dict, annotator=None, gene_runner=None) -> Path:
    """Annotate once, run gene analysis per group; return the manifest path."""
    paths = config["paths"]
    metal_dir = Path(paths["metal_dir"])
    out_dir = Path(paths["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    magma_bin = Path(config["magma_bin"])
    write_run_provenance(
        out_dir, tool="MAGMA", binary=magma_bin,
        version=magma_runner_mod.tool_version(magma_bin),
        config_path=config.get("_config_path"), repo_root=REPO_ROOT,
        extra={"bfile": config.get("bfile"), "gene_loc": config.get("gene_loc")},
    )

    if annotator is None:
        annotator = _default_annotator(config, out_dir)
    annot = annotator()
    if gene_runner is None:
        gene_runner = _default_gene_runner(config, annot)

    gene_frames: list[pd.DataFrame] = []
    sig_frames: list[pd.DataFrame] = []
    manifest_rows: list[dict] = []

    for group in config["groups"]:
        tbl = metal_dir / f"{group}_1.tbl"
        if not tbl.exists():
            manifest_rows.append({"group": group, "n_snps_in": 0, "n_genes_tested": 0,
                                  "n_genes_sig": 0, "status": "missing_input"})
            logger.warning("[%s] missing %s", group, tbl)
            continue
        pval_file = out_dir / f"{group}.pval.txt"
        n_snps = prep_pval_input(tbl, pval_file)
        result = gene_runner(group, pval_file, out_dir)
        if result is None or len(result) == 0:
            manifest_rows.append({"group": group, "n_snps_in": n_snps,
                                  "n_genes_tested": 0, "n_genes_sig": 0,
                                  "status": "no_genes"})
            continue
        gene_frames.append(result)
        # Bonferroni across genes tested in this group.
        threshold = 0.05 / len(result)
        sig = result[pd.to_numeric(result["pval"], errors="coerce") < threshold]
        sig_frames.append(sig)
        manifest_rows.append({"group": group, "n_snps_in": n_snps,
                              "n_genes_tested": len(result), "n_genes_sig": len(sig),
                              "status": "analyzed"})
        logger.info("[%s] %d SNPs -> %d genes, %d significant",
                    group, n_snps, len(result), len(sig))

    (pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()) \
        .to_csv(out_dir / "gene_results.csv", index=False)
    (pd.concat(sig_frames, ignore_index=True) if sig_frames else pd.DataFrame()) \
        .to_csv(out_dir / "significant_genes.csv", index=False)
    manifest_path = out_dir / "magma_manifest.csv"
    pd.DataFrame(manifest_rows, columns=_MANIFEST_COLUMNS).to_csv(manifest_path, index=False)
    logger.info("Wrote gene results, significant genes, manifest -> %s", out_dir)
    return manifest_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = p.parse_args()
    with open(args.config) as f:
        config = yaml.safe_load(f)
    run(config)


if __name__ == "__main__":
    main()
