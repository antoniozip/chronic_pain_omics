from pathlib import Path

import pandas as pd

from cp_multiomics.genomics import magma_runner
from cp_multiomics.genomics.magma_runner import (
    parse_genes_out,
    prep_pval_input,
    write_snploc,
)

_FIXTURE = Path(__file__).parent / "fixtures" / "magma_sample.genes.out"
_TBL_HEADER = ("MarkerName\tAllele1\tAllele2\tWeight\tZscore\tP-value\tDirection\t"
               "HetISq\tHetChiSq\tHetDf\tHetPVal\n")


def test_write_snploc(tmp_path):
    bim = tmp_path / "ref.bim"
    bim.write_text("1\trs1\t0\t10539\tA\tC\n2\trs2\t0\t20000\tG\tT\n")
    out = tmp_path / "snploc.txt"
    result = write_snploc(bim, out)
    assert result == out
    df = pd.read_csv(out, sep="\t", header=None, names=["SNP", "CHR", "BP"])
    assert df["SNP"].tolist() == ["rs1", "rs2"]
    assert df["CHR"].tolist() == [1, 2]
    assert df["BP"].tolist() == [10539, 20000]
    # no header line
    assert out.read_text().splitlines()[0].split("\t")[0] == "rs1"


def test_prep_pval_input_drops_bad_rows(tmp_path):
    tbl = tmp_path / "grp_1.tbl"
    tbl.write_text(
        _TBL_HEADER
        + "rs1\ta\tg\t715631.00\t5.5\t1.77e-08\t++++\t0.0\t1\t3\t0.9\n"
        + "rs2\ta\tg\tNA\t0.1\tNA\t+?-?\t0.0\t0\t0\t1\n"       # NaN P and N -> drop
    )
    out = tmp_path / "pval.txt"
    n = prep_pval_input(tbl, out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["SNP", "P", "N"]
    assert df["SNP"].tolist() == ["rs1"]
    assert df["N"].tolist() == [715631]


def test_parse_genes_out():
    df = parse_genes_out(_FIXTURE, "back_pain")
    assert list(df.columns) == [
        "group", "gene", "chr", "start", "stop", "n_snps",
        "n_param", "n", "zstat", "pval",
    ]
    assert set(df["group"]) == {"back_pain"}
    assert len(df) == 3
    row = df.set_index("gene").loc["79501"]
    assert row["chr"] == 1
    assert row["n_snps"] == 12
    assert row["zstat"] == 5.512
    assert row["pval"] == 1.77e-08


def test_annotate_appends_suffix_to_a_dotted_prefix(tmp_path, monkeypatch):
    """Path.with_suffix would REPLACE '.0', pointing at a file MAGMA never wrote."""
    monkeypatch.setattr(magma_runner, "_run", lambda cmd: None)
    prefix = tmp_path / "chronic_pain_v1.0"
    got = magma_runner.annotate(Path("magma"), tmp_path / "s.loc",
                                tmp_path / "g.loc", prefix)
    assert got.name == "chronic_pain_v1.0.genes.annot"


def test_gene_analysis_appends_suffix_to_a_dotted_prefix(tmp_path, monkeypatch):
    monkeypatch.setattr(magma_runner, "_run", lambda cmd: None)
    prefix = tmp_path / "chronic_pain_v1.0"
    got = magma_runner.run_gene_analysis(Path("magma"), tmp_path / "ref",
                                         tmp_path / "p.txt",
                                         tmp_path / "a.annot", prefix)
    assert got.name == "chronic_pain_v1.0.genes.out"


def test_undotted_prefixes_are_unchanged(tmp_path, monkeypatch):
    """The production group slugs have no dot; behaviour must not shift."""
    monkeypatch.setattr(magma_runner, "_run", lambda cmd: None)
    prefix = tmp_path / "back_pain"
    assert magma_runner.annotate(Path("m"), tmp_path / "s", tmp_path / "g",
                                 prefix).name == "back_pain.genes.annot"
    assert magma_runner.run_gene_analysis(Path("m"), tmp_path / "r", tmp_path / "p",
                                          tmp_path / "a",
                                          prefix).name == "back_pain.genes.out"


def test_tool_version_captures_the_banner(tmp_path):
    stub = tmp_path / "magma"
    stub.write_text("#!/bin/sh\necho 'MAGMA version: v1.10 (linux)'\n")
    stub.chmod(0o755)
    assert "v1.10" in magma_runner.tool_version(stub)


def test_tool_version_returns_unknown_rather_than_raising(tmp_path):
    stub = tmp_path / "magma"
    stub.write_text("#!/bin/sh\nexit 2\n")
    stub.chmod(0o755)
    assert magma_runner.tool_version(stub) == "unknown"


def test_tool_version_returns_unknown_for_a_missing_binary(tmp_path):
    assert magma_runner.tool_version(tmp_path / "nope") == "unknown"
