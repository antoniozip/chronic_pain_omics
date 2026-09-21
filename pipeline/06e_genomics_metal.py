"""Pool the harmonized chronic-pain GWAS into cross-cohort METAL meta-analyses.

Consumes the Stage-1 harmonized per-study tables + manifest, assigns each fetched
study to a broad clinical pain group, collapses UK Biobank slices to one
representative per cohort (cohort_expansion_plan.md §7.7), and runs the real METAL
binary once per poolable group (>= min_cohorts independent cohorts). Single-cohort
groups are reported without pooling; non-pain phenotypes are excluded. Every study
is accounted for in a PRISMA-style grouping manifest.

    uv run python pipeline/06e_genomics_metal.py

Retires the manuscript's "precluding standard METAL-based meta-analysis" claim.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

from cp_multiomics.genomics import metal_runner as metal_runner_mod
from cp_multiomics.genomics.grouping import (
    poolable_groups,
    select_representatives,
)
from cp_multiomics.genomics.metal_runner import (
    parse_metal_tbl,
    run_metal,
    write_metal_script,
)
from cp_multiomics.genomics.provenance import write_run_provenance

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "conf" / "genomics" / "metal.yaml"

_GROUPING_COLUMNS = [
    "accession", "phenotype", "group", "cohort_family",
    "is_representative", "pooled_status", "n",
]


def _default_metal_runner(metal_bin: Path):
    """Build the production runner: write script -> run METAL -> parse .tbl."""
    def runner(group: str, table_paths: list[Path], out_dir: Path) -> pd.DataFrame | None:
        out_prefix = out_dir / f"{group}_"
        script_path = out_dir / f"{group}.metal"
        write_metal_script(group, table_paths, out_prefix, script_path)
        run_metal(script_path, metal_bin)
        tbl_path = out_dir / f"{group}_1.tbl"   # METAL appends '1' to the prefix
        if not tbl_path.exists():
            logger.error("[%s] METAL produced no output at %s", group, tbl_path)
            return None
        return parse_metal_tbl(tbl_path, group)
    return runner


def run(config: dict, metal_runner=None) -> Path:
    """Assign groups, collapse per §7.7, pool poolable groups; return manifest path."""
    paths = config["paths"]
    interim_dir = Path(paths["interim_dir"])
    out_dir = Path(paths["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    min_cohorts = int(config.get("min_cohorts", 2))
    gws_threshold = float(config.get("gws_threshold", 5e-8))
    if metal_runner is None:
        metal_runner = _default_metal_runner(Path(config["metal_bin"]))

    metal_bin = Path(config["metal_bin"])
    write_run_provenance(
        out_dir, tool="METAL", binary=metal_bin,
        version=metal_runner_mod.tool_version(metal_bin),
        config_path=config.get("_config_path"), repo_root=REPO_ROOT,
        extra={"min_cohorts": min_cohorts, "gws_threshold": gws_threshold},
    )

    manifest = pd.read_csv(paths["sumstats_manifest"])
    fetched = manifest[manifest["disposition"] == "fetched"].copy()
    logger.info("Grouping %d fetched studies", len(fetched))

    rep = select_representatives(fetched)
    grouping_path = out_dir / "grouping_manifest.csv"
    rep[_GROUPING_COLUMNS].to_csv(grouping_path, index=False)

    counts = poolable_groups(rep)
    pooled_frames: list[pd.DataFrame] = []
    status_rows: list[dict] = []

    for group, n_cohorts in sorted(counts.items()):
        reps = rep[(rep["group"] == group) & rep["is_representative"]]
        if n_cohorts < min_cohorts:
            status_rows.append({"group": group, "n_cohorts": n_cohorts,
                                "status": "not_pooled_single_cohort",
                                "n_markers": 0, "n_gws": 0})
            logger.info("[%s] single cohort -> not pooled", group)
            continue
        table_paths = [interim_dir / f"{acc}.tsv" for acc in reps["accession"]]
        # METAL warns and exits 0 on a missing PROCESS file, so a group would
        # pool fewer cohorts than group_status.csv reports, with nothing in the
        # logs to say so. Refuse the group instead.
        absent = [p for p in table_paths if not p.exists()]
        if absent:
            raise FileNotFoundError(
                f"[{group}] missing harmonized table(s): "
                f"{', '.join(p.name for p in absent)}"
            )
        result = metal_runner(group, table_paths, out_dir)
        n_markers = n_gws = 0
        if result is not None:
            n_markers = len(result)
            # Full per-marker results live in the per-group .tbl; the summary
            # keeps only genome-wide-significant hits so it stays usable.
            gws = result[pd.to_numeric(result["pval"], errors="coerce") < gws_threshold]
            n_gws = len(gws)
            pooled_frames.append(gws)
        status_rows.append({"group": group, "n_cohorts": n_cohorts,
                            "status": "pooled", "n_markers": n_markers, "n_gws": n_gws})
        logger.info("[%s] pooled %d cohorts -> %d markers, %d genome-wide-sig",
                    group, n_cohorts, n_markers, n_gws)

    # An empty CSV with no header cannot be read back by any consumer, so seed
    # the empty case with the METAL output schema.
    pooled = (pd.concat(pooled_frames, ignore_index=True)
              if pooled_frames
              else pd.DataFrame(columns=metal_runner_mod.OUTPUT_COLUMNS))
    pooled.to_csv(out_dir / "pooled_summary.csv", index=False)
    pd.DataFrame(status_rows).to_csv(out_dir / "group_status.csv", index=False)
    logger.info("Wrote grouping manifest, pooled summary, group status -> %s", out_dir)
    return grouping_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = p.parse_args()
    with open(args.config) as f:
        config = yaml.safe_load(f)
    run(config)


if __name__ == "__main__":
    main()
