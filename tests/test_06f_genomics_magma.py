import csv
import importlib.util
from pathlib import Path

import pandas as pd

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_06f", Path(__file__).parent.parent / "pipeline" / "06f_genomics_magma.py"
)
pipeline_06f = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_06f)

_TBL_HEADER = ("MarkerName\tAllele1\tAllele2\tWeight\tZscore\tP-value\tDirection\t"
               "HetISq\tHetChiSq\tHetDf\tHetPVal\n")


def _cfg(tmp_path):
    metal_dir = tmp_path / "metal"
    metal_dir.mkdir()
    for grp in ("back_pain", "fibromyalgia"):
        (metal_dir / f"{grp}_1.tbl").write_text(
            _TBL_HEADER
            + "rs1\ta\tg\t700000\t5.5\t1.77e-08\t++++\t0.0\t1\t3\t0.9\n"
            + "rs2\tc\tt\t700000\t0.3\t0.7\t+-+-\t0.0\t1\t3\t0.9\n"
        )
    return {
        "paths": {"metal_dir": str(metal_dir), "out_dir": str(tmp_path / "out")},
        "magma_bin": "unused", "bfile": "unused", "gene_loc": "unused",
        "groups": ["back_pain", "fibromyalgia"],
    }


def test_orchestrator_runs_each_group_and_flags_significant(tmp_path):
    calls = []

    def stub_annotator():
        return Path("dummy.genes.annot")

    def stub_gene_runner(group, pval_file, out_dir):
        calls.append(group)
        return pd.DataFrame([
            {"group": group, "gene": "G_SIG", "chr": 1, "start": 1, "stop": 2,
             "n_snps": 10, "n_param": 3, "n": 700000, "zstat": 5.6, "pval": 1e-8},
            {"group": group, "gene": "G_NULL1", "chr": 1, "start": 3, "stop": 4,
             "n_snps": 8, "n_param": 3, "n": 700000, "zstat": 1.2, "pval": 0.11},
            {"group": group, "gene": "G_NULL2", "chr": 2, "start": 5, "stop": 6,
             "n_snps": 6, "n_param": 3, "n": 700000, "zstat": -0.3, "pval": 0.63},
        ])

    cfg = _cfg(tmp_path)
    manifest_path = pipeline_06f.run(cfg, annotator=stub_annotator,
                                     gene_runner=stub_gene_runner)

    assert calls == ["back_pain", "fibromyalgia"]        # once per group

    out_dir = Path(cfg["paths"]["out_dir"])
    gene_results = pd.read_csv(out_dir / "gene_results.csv")
    assert len(gene_results) == 6                        # 3 genes x 2 groups

    sig = pd.read_csv(out_dir / "significant_genes.csv")
    # Bonferroni 0.05/3 = 0.0167 -> only G_SIG passes, in both groups
    assert set(sig["gene"]) == {"G_SIG"}
    assert len(sig) == 2

    manifest = {r["group"]: r for r in csv.DictReader(open(manifest_path))}
    assert set(manifest) == {"back_pain", "fibromyalgia"}
    assert manifest["back_pain"]["n_genes_tested"] == "3"
    assert manifest["back_pain"]["n_genes_sig"] == "1"
    assert manifest["back_pain"]["n_snps_in"] == "2"
    assert manifest["back_pain"]["status"] == "analyzed"
