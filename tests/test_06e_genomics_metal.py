import csv
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from cp_multiomics.genomics import metal_runner

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_06e", Path(__file__).parent.parent / "pipeline" / "06e_genomics_metal.py"
)
pipeline_06e = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_06e)


def _cfg(tmp_path, min_cohorts=2):
    manifest = tmp_path / "sumstats_manifest.csv"
    rows = [
        # back_pain: 1 UKB + 1 independent -> 2 cohorts -> pooled
        {"accession": "GCST_BP1", "phenotype": "back pain", "cohort": "UKB",
         "n": 300000, "disposition": "fetched"},
        {"accession": "GCST_BP2", "phenotype": "chronic back pain", "cohort": "",
         "n": 50000, "disposition": "fetched"},
        # fibromyalgia: 2 UKB slices -> collapse to 1 cohort -> single-cohort
        {"accession": "GCST_FM1", "phenotype": "fibromyalgia", "cohort": "UKB",
         "n": 400000, "disposition": "fetched"},
        {"accession": "GCST_FM2", "phenotype": "icd10 m79: fibromyalgia",
         "cohort": "PRESUMED_BIOBANK", "n": 380000, "disposition": "fetched"},
        # nonpain
        {"accession": "GCST_NP1", "phenotype": "menstruation quality of life impact",
         "cohort": "", "n": 5734, "disposition": "fetched"},
        # excluded upstream -> must be ignored
        {"accession": "GCST_X", "phenotype": "back pain", "cohort": "UKB",
         "n": 1, "disposition": "excluded"},
    ]
    with open(manifest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    # Every fetched study must have a harmonized table on disk; the orchestrator
    # refuses to pool a group whose inputs are missing.
    interim = tmp_path / "interim"
    interim.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if row["disposition"] == "fetched":
            (interim / f"{row['accession']}.tsv").write_text("rsid\tbeta\n")
    return {
        "paths": {
            "sumstats_manifest": str(manifest),
            "interim_dir": str(interim),
            "out_dir": str(tmp_path / "out"),
        },
        "metal_bin": "unused-in-test",
        "min_cohorts": min_cohorts,
    }


def test_orchestrator_pools_only_multicohort_groups(tmp_path):
    calls = []

    def stub_runner(group, table_paths, out_dir):
        calls.append((group, [Path(p).name for p in table_paths]))
        return pd.DataFrame([
            {"group": group, "rsid": "rs_gws", "allele1": "a", "allele2": "g",
             "weight": 350000.0, "zscore": 6.0, "pval": 1e-9, "direction": "++",
             "hetisq": 0.0, "hetpval": 0.9, "n_studies": 2},
            {"group": group, "rsid": "rs_null", "allele1": "c", "allele2": "t",
             "weight": 350000.0, "zscore": 0.5, "pval": 0.6, "direction": "+-",
             "hetisq": 0.0, "hetpval": 0.9, "n_studies": 2},
        ])

    cfg = _cfg(tmp_path)
    manifest_path = pipeline_06e.run(cfg, metal_runner=stub_runner)

    # grouping manifest: every fetched study present with a pooled_status
    gm = {r["accession"]: r for r in csv.DictReader(open(manifest_path))}
    assert set(gm) == {"GCST_BP1", "GCST_BP2", "GCST_FM1", "GCST_FM2", "GCST_NP1"}
    assert gm["GCST_BP1"]["group"] == "back_pain"
    assert gm["GCST_NP1"]["pooled_status"] == "nonpain"
    # fibromyalgia UKB slices: larger n representative, other collapsed
    assert gm["GCST_FM1"]["is_representative"] == "True"
    assert gm["GCST_FM2"]["pooled_status"] == "collapsed_overlap"

    # METAL invoked for back_pain (2 cohorts) only
    pooled_groups = [c[0] for c in calls]
    assert pooled_groups == ["back_pain"]
    # both back_pain representatives passed as tables
    assert set(calls[0][1]) == {"GCST_BP1.tsv", "GCST_BP2.tsv"}

    # pooled_summary holds ONLY genome-wide-significant markers (p<5e-8)
    ps = pd.read_csv(Path(cfg["paths"]["out_dir"]) / "pooled_summary.csv")
    assert ps["rsid"].tolist() == ["rs_gws"]          # rs_null (p=0.6) filtered out
    assert (ps["group"] == "back_pain").all()

    # group_status: fibromyalgia single-cohort; back_pain totals + gws count
    gs = {r["group"]: r for r in csv.DictReader(
        open(Path(cfg["paths"]["out_dir"]) / "group_status.csv"))}
    assert gs["back_pain"]["status"] == "pooled"
    assert gs["back_pain"]["n_markers"] == "2"
    assert gs["back_pain"]["n_gws"] == "1"
    assert gs["fibromyalgia"]["status"] == "not_pooled_single_cohort"


def test_run_raises_when_a_representative_table_is_missing(tmp_path):
    """METAL warns and exits 0 on a missing PROCESS file, silently pooling
    fewer cohorts than group_status.csv reports. Fail loudly instead."""
    cfg = _cfg(tmp_path)
    (Path(cfg["paths"]["interim_dir"]) / "GCST_BP2.tsv").unlink()

    def stub_runner(group, table_paths, out_dir):  # pragma: no cover - must not run
        raise AssertionError("METAL must not be invoked with a missing input")

    with pytest.raises(FileNotFoundError, match="GCST_BP2"):
        pipeline_06e.run(cfg, metal_runner=stub_runner)


def test_pooled_summary_has_a_header_when_no_group_pools(tmp_path):
    """An empty CSV with no header is unreadable by every downstream consumer."""
    cfg = _cfg(tmp_path, min_cohorts=99)   # nothing reaches the cohort threshold

    def stub_runner(group, table_paths, out_dir):  # pragma: no cover - must not run
        raise AssertionError("no group should pool")

    pipeline_06e.run(cfg, metal_runner=stub_runner)
    out = pd.read_csv(Path(cfg["paths"]["out_dir"]) / "pooled_summary.csv")
    assert out.empty
    assert list(out.columns) == metal_runner.OUTPUT_COLUMNS
