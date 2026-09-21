from pathlib import Path

from cp_multiomics.genomics import metal_runner
from cp_multiomics.genomics.metal_runner import parse_metal_tbl, write_metal_script

_FIXTURE = Path(__file__).parent / "fixtures" / "metal_sample.tbl"


def test_write_metal_script_structure(tmp_path):
    tables = [tmp_path / f"s{i}.tsv" for i in range(3)]
    for t in tables:
        t.write_text("rsid\teffect_allele\tother_allele\tbeta\tpval\tn\n")
    out_prefix = tmp_path / "back_pain_"
    script_path = tmp_path / "back_pain.metal"

    result = write_metal_script("back_pain", tables, out_prefix, script_path)
    assert result == script_path
    text = script_path.read_text()

    assert "SCHEME SAMPLESIZE" in text
    assert "SEPARATOR TAB" in text
    assert "MARKER rsid" in text
    assert "ALLELE effect_allele other_allele" in text
    assert "EFFECT beta" in text
    assert "PVALUE pval" in text
    assert "WEIGHT n" in text
    process_lines = [ln for ln in text.splitlines() if ln.startswith("PROCESS ")]
    assert len(process_lines) == 3
    assert [ln.split(" ", 1)[1] for ln in process_lines] == [str(t) for t in tables]
    assert f"OUTFILE {out_prefix} .tbl" in text
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    assert lines[-2] == "ANALYZE HETEROGENEITY"
    assert lines[-1] == "QUIT"


def test_parse_metal_tbl():
    df = parse_metal_tbl(_FIXTURE, "back_pain")
    assert list(df.columns) == [
        "group", "rsid", "allele1", "allele2", "weight", "zscore",
        "pval", "direction", "hetisq", "hetpval", "n_studies",
    ]
    assert set(df["group"]) == {"back_pain"}
    assert len(df) == 2
    row = df.set_index("rsid").loc["rs1"]
    assert row["zscore"] == 4.324
    assert row["pval"] == 1.532e-05
    assert row["direction"] == "++"
    assert row["n_studies"] == 2       # two '+' contributions, no '?'


def test_parse_metal_tbl_counts_present_studies(tmp_path):
    tbl = tmp_path / "x.tbl"
    tbl.write_text(
        "MarkerName\tAllele1\tAllele2\tWeight\tZscore\tP-value\tDirection\t"
        "HetISq\tHetChiSq\tHetDf\tHetPVal\n"
        "rs9\ta\tg\t9000.00\t2.10\t0.036\t+-?\t0.0\t0.10\t2\t0.95\n"
    )
    df = parse_metal_tbl(tbl, "headache")
    assert df.iloc[0]["n_studies"] == 2   # '+' and '-' count, '?' does not


def test_tool_version_captures_the_banner(tmp_path):
    stub = tmp_path / "metal"
    stub.write_text("#!/bin/sh\necho 'MetaAnalysis Helper - (c) 2007-2009 Goncalo Abecasis'\n"
                    "echo 'This version released on 2011-03-25'\n")
    stub.chmod(0o755)
    assert "MetaAnalysis Helper" in metal_runner.tool_version(stub)


def test_tool_version_returns_unknown_rather_than_raising(tmp_path):
    """A version probe must never break an analysis."""
    stub = tmp_path / "metal"
    stub.write_text("#!/bin/sh\nexit 1\n")
    stub.chmod(0o755)
    assert metal_runner.tool_version(stub) == "unknown"


def test_tool_version_returns_unknown_for_a_missing_binary(tmp_path):
    assert metal_runner.tool_version(tmp_path / "does-not-exist") == "unknown"
