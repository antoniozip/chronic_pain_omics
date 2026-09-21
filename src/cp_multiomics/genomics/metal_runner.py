"""Write METAL control scripts, run the METAL binary, and parse its pooled output.

The pool uses METAL's SAMPLESIZE (weighted-sum-of-Z) scheme because the cohort
mixes linear-model betas and case/control log-ORs, which share no common
effect-size scale: SAMPLESIZE combines them from effect direction, p-value, and N.
Effect sizes still enter as the EFFECT column so METAL derives each study's sign
(and flips it consistently when it normalizes alleles to a reference orientation).
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pandas as pd

from cp_multiomics.genomics.provenance import probe_version

logger = logging.getLogger(__name__)

# METAL SAMPLESIZE + ANALYZE HETEROGENEITY output columns -> our tidy names.
_TBL_RENAME = {
    "MarkerName": "rsid",
    "Allele1": "allele1",
    "Allele2": "allele2",
    "Weight": "weight",
    "Zscore": "zscore",
    "P-value": "pval",
    "Direction": "direction",
    "HetISq": "hetisq",
    "HetPVal": "hetpval",
}
OUTPUT_COLUMNS = [
    "group", "rsid", "allele1", "allele2", "weight", "zscore",
    "pval", "direction", "hetisq", "hetpval", "n_studies",
]


def tool_version(metal_bin: Path) -> str:
    """Return METAL's startup banner, or "unknown".

    METAL has no --version flag: it prints its banner and then reads a command
    script from stdin, so the probe feeds it EOF immediately.
    """
    return probe_version(metal_bin)


def write_metal_script(
    group: str,
    table_paths: list[Path],
    out_prefix: Path,
    script_path: Path,
) -> Path:
    """Write a SAMPLESIZE-scheme METAL control script; return its path."""
    lines = [
        f"# METAL meta-analysis for pain group: {group}",
        "SCHEME SAMPLESIZE",
        "SEPARATOR TAB",
        "MARKER rsid",
        "ALLELE effect_allele other_allele",
        "EFFECT beta",
        "PVALUE pval",
        "WEIGHT n",
    ]
    lines += [f"PROCESS {path}" for path in table_paths]
    lines += [f"OUTFILE {out_prefix} .tbl", "ANALYZE HETEROGENEITY", "QUIT", ""]
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(lines))
    return script_path


def run_metal(script_path: Path, metal_bin: Path) -> None:
    """Invoke the METAL binary on a control script; raise on failure."""
    result = subprocess.run(
        [str(metal_bin), str(script_path)],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"METAL failed ({result.returncode}) on {script_path}:\n{result.stderr}"
        )
    logger.info("METAL completed for %s", script_path.name)


def parse_metal_tbl(tbl_path: Path, group: str) -> pd.DataFrame:
    """Parse a METAL SAMPLESIZE+heterogeneity .tbl into a tidy per-marker frame."""
    df = pd.read_csv(tbl_path, sep="\t")
    df = df.rename(columns=_TBL_RENAME)
    df["group"] = group
    # n_studies = studies actually contributing a marker: '+'/'-' in Direction,
    # excluding '?' (marker absent from that study).
    df["n_studies"] = df["direction"].astype(str).apply(
        lambda d: sum(c in "+-" for c in d)
    )
    for col in ("weight", "zscore", "pval", "hetisq", "hetpval"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df[OUTPUT_COLUMNS]
