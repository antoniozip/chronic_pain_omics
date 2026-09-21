"""Prepare inputs for and parse outputs from MAGMA gene-based analysis.

MAGMA aggregates SNP-level association signal to genes, accounting for LD via a
reference panel. It complements the SNP-level METAL meta-analysis: genes can reach
significance from many sub-threshold SNPs that no single marker would surface.

SNPs are matched between our (GRCh38) summary statistics and the (GRCh37) 1000G
reference by rsID, which is build-stable, so the coordinate build differs but the
analysis stays consistent on the reference side. Per-SNP N is used (meta-analysis
sample size varies by SNP) via MAGMA's ncol option.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pandas as pd

from cp_multiomics.genomics.provenance import probe_version

logger = logging.getLogger(__name__)

_GENES_OUT_RENAME = {
    "GENE": "gene", "CHR": "chr", "START": "start", "STOP": "stop",
    "NSNPS": "n_snps", "NPARAM": "n_param", "N": "n", "ZSTAT": "zstat", "P": "pval",
}
_GENES_OUT_COLUMNS = [
    "group", "gene", "chr", "start", "stop", "n_snps",
    "n_param", "n", "zstat", "pval",
]


def tool_version(magma_bin: Path) -> str:
    """Return MAGMA's --version banner, or "unknown"."""
    return probe_version(magma_bin, args=["--version"])


def _with_added_suffix(prefix: Path, suffix: str) -> Path:
    """Append `suffix` to a path used as a MAGMA `--out` prefix.

    Path.with_suffix would *replace* an existing extension, so a prefix
    containing a dot (e.g. "chronic_pain_v1.0") would resolve to a different
    file than the one MAGMA actually writes.
    """
    return prefix.with_name(prefix.name + suffix)


def write_snploc(bim_path: Path, out_path: Path) -> Path:
    """Write a MAGMA SNP-location file (SNP CHR BP, no header) from a plink .bim."""
    bim = pd.read_csv(
        bim_path, sep="\t", header=None,
        names=["chr", "snp", "cm", "bp", "a1", "a2"],
        usecols=["chr", "snp", "bp"], dtype={"snp": str},
    )
    bim[["snp", "chr", "bp"]].to_csv(out_path, sep="\t", header=False, index=False)
    return out_path


def prep_pval_input(tbl_path: Path, out_path: Path) -> int:
    """Write a MAGMA p-value input (SNP P N) from a METAL .tbl; return rows written."""
    tbl = pd.read_csv(tbl_path, sep="\t", usecols=["MarkerName", "P-value", "Weight"])
    out = pd.DataFrame({
        "SNP": tbl["MarkerName"].astype(str),
        "P": pd.to_numeric(tbl["P-value"], errors="coerce"),
        "N": pd.to_numeric(tbl["Weight"], errors="coerce"),
    })
    out = out.dropna(subset=["P", "N"])
    out["N"] = out["N"].round().astype(int)
    out.to_csv(out_path, sep="\t", index=False)
    return len(out)


def annotate(magma_bin: Path, snploc: Path, gene_loc: Path, out_prefix: Path) -> Path:
    """Run MAGMA SNP->gene annotation; return the .genes.annot path."""
    _run([str(magma_bin), "--annotate", "--snp-loc", str(snploc),
          "--gene-loc", str(gene_loc), "--out", str(out_prefix)])
    return _with_added_suffix(out_prefix, ".genes.annot")


def run_gene_analysis(
    magma_bin: Path, bfile: Path, pval_file: Path, gene_annot: Path,
    out_prefix: Path, n_col: str = "N",
) -> Path:
    """Run MAGMA gene analysis; return the .genes.out path."""
    _run([str(magma_bin), "--bfile", str(bfile),
          "--pval", str(pval_file), f"ncol={n_col}",
          "--gene-annot", str(gene_annot), "--out", str(out_prefix)])
    return _with_added_suffix(out_prefix, ".genes.out")


def parse_genes_out(path: Path, group: str) -> pd.DataFrame:
    """Parse a MAGMA .genes.out (whitespace-delimited) into a tidy frame."""
    df = pd.read_csv(path, sep=r"\s+", dtype={"GENE": str})
    df = df.rename(columns=_GENES_OUT_RENAME)
    df["group"] = group
    return df[_GENES_OUT_COLUMNS]


def _run(cmd: list[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"MAGMA failed ({result.returncode}): {' '.join(cmd)}\n"
                           f"{result.stdout[-2000:]}\n{result.stderr[-2000:]}")
    logger.info("MAGMA ok: %s", cmd[1] if len(cmd) > 1 else cmd[0])
