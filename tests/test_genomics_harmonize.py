import gzip
import math
from pathlib import Path

import pandas as pd
import pytest

from cp_multiomics.genomics.harmonize import _resolve_rsid, harmonize_sumstats
from cp_multiomics.genomics.resolve import SumstatsStudy

_HEADER = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
          "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t" \
          "hm_coordinate_conversion\thm_code\trsid\n"

_OR_HEADER = "chromosome\tbase_pair_location\teffect_allele\tother_allele\todds_ratio\t" \
             "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t" \
             "hm_coordinate_conversion\thm_code\trsid\tci_upper\tci_lower\n"

_NO_EFFECT_HEADER = "chromosome\tbase_pair_location\teffect_allele\tother_allele\t" \
                     "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t" \
                     "hm_coordinate_conversion\thm_code\trsid\n"

_QNORM_975 = 1.959964


def _write_gz(tmp_path: Path, header: str, rows: list[str], name: str = "in.h.tsv.gz") -> Path:
    p = tmp_path / name
    with gzip.open(p, "wt") as f:
        f.write(header)
        for r in rows:
            f.write(r + "\n")
    return p


def _study():
    return SumstatsStudy(accession="GCST90091914", trait="Fibromyalgia",
                         sample_size_text="", n=45409, cohort="",
                         is_burden=False, phenotype="fibromyalgia")


def test_columns_are_remapped_and_n_attached(tmp_path):
    gz = _write_gz(tmp_path, _HEADER, [
        "1\t10177\tAC\tA\t-0.0877\t0.2043\t0.4573\t0.6678\tNA\tlo\t11\trs367896724",
    ])
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["rsid", "chromosome", "position", "effect_allele",
                                "other_allele", "beta", "se", "eaf", "pval", "n",
                                "study_id", "cohort", "phenotype", "build", "effect_unit"]
    row = df.iloc[0]
    assert row["rsid"] == "rs367896724"
    assert row["position"] == 10177
    assert row["se"] == 0.2043
    assert row["n"] == 45409
    assert row["study_id"] == "GCST90091914"
    assert row["phenotype"] == "fibromyalgia"
    assert row["effect_unit"] == "beta"


def test_rows_without_rsid_or_se_are_dropped(tmp_path):
    gz = _write_gz(tmp_path, _HEADER, [
        "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1",   # keep
        "1\t10178\tG\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\tNA",     # no rsid -> drop
        "1\t10179\tG\tA\tNA\tNA\t0.45\t0.66\tNA\tlo\t11\trs3",         # no beta/se -> drop
    ])
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert df["rsid"].tolist() == ["rs1"]


def test_or_only_row_derives_beta_and_se_from_ci(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t2.0\tNA\t0.45\t0.66\tNA\tlo\t11\trs1\t3.0\t1.333333",
    ], name="or.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    row = df.iloc[0]
    assert row["beta"] == pytest.approx(math.log(2.0))
    expected_se = (math.log(3.0) - math.log(1.333333)) / (2 * _QNORM_975)
    assert row["se"] == pytest.approx(expected_se)
    assert row["effect_unit"] == "log_OR"


def test_or_only_row_with_populated_se_uses_se_not_ci(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t2.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs1\t3.0\t1.333333",
    ], name="or_se.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    row = df.iloc[0]
    assert row["se"] == pytest.approx(0.15)
    assert row["effect_unit"] == "log_OR"


def test_no_beta_no_odds_ratio_raises(tmp_path):
    gz = _write_gz(tmp_path, _NO_EFFECT_HEADER, [
        "1\t10177\tAC\tA\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1",
    ], name="no_effect.h.tsv.gz")
    out = tmp_path / "out.tsv"
    with pytest.raises(ValueError):
        harmonize_sumstats(gz, _study(), out)


def test_or_only_row_without_se_or_ci_is_dropped(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t2.0\tNA\t0.45\t0.66\tNA\tlo\t11\trs1\tNA\tNA",
    ], name="or_no_ci.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 0


def test_odds_ratio_zero_is_dropped(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t0.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs1\tNA\tNA",
        "1\t10178\tAC\tA\t2.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs2\tNA\tNA",
    ], name="or_zero.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert df["rsid"].tolist() == ["rs2"]
    import numpy as np
    assert np.isinf(df["beta"]).sum() == 0


def test_odds_ratio_negative_is_dropped(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t-1.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs1\tNA\tNA",
        "1\t10178\tAC\tA\t2.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs2\tNA\tNA",
    ], name="or_negative.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert df["rsid"].tolist() == ["rs2"]


def test_all_nan_beta_column_falls_through_to_odds_ratio(tmp_path):
    header = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
             "odds_ratio\tstandard_error\teffect_allele_frequency\tp_value\t" \
             "variant_id\thm_coordinate_conversion\thm_code\trsid\n"
    gz = _write_gz(tmp_path, header, [
        "1\t10177\tAC\tA\tNA\t2.0\t0.15\t0.45\t0.66\tNA\tlo\t11\trs1",
    ], name="beta_all_nan.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    row = df.iloc[0]
    assert row["beta"] == pytest.approx(math.log(2.0))
    assert row["effect_unit"] == "log_OR"


def test_hm_rsid_used_when_bare_rsid_column_absent(tmp_path):
    # Fully-harmonised GWAS-SSF files name the identifier `hm_rsid` and carry no
    # bare `rsid` column; 31 of the 61 kept chronic-pain studies are this shape.
    header = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
             "standard_error\teffect_allele_frequency\tp_value\tvariant_id\thm_rsid\n"
    gz = _write_gz(tmp_path, header, [
        "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\t1_10177_AC_A\trs999",
    ], name="hm_rsid.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["rsid"] == "rs999"


def test_variant_id_used_when_no_rsid_columns(tmp_path):
    header = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
             "standard_error\teffect_allele_frequency\tp_value\tvariant_id\n"
    gz = _write_gz(tmp_path, header, [
        "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\t1_10177_AC_A",
    ], name="variant_only.h.tsv.gz")
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 1
    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["rsid"] == "1_10177_AC_A"


def test_bare_rsid_preferred_over_hm_rsid(tmp_path):
    # Files that already harmonized cleanly used the bare `rsid` column; keep
    # that stable when both are present so existing outputs do not change.
    header = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
             "standard_error\teffect_allele_frequency\tp_value\trsid\thm_rsid\n"
    gz = _write_gz(tmp_path, header, [
        "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\trs_bare\trs_hm",
    ], name="both_rsid.h.tsv.gz")
    out = tmp_path / "out.tsv"
    harmonize_sumstats(gz, _study(), out)
    df = pd.read_csv(out, sep="\t")
    assert df.iloc[0]["rsid"] == "rs_bare"


def test_no_identifier_column_raises(tmp_path):
    header = "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t" \
             "standard_error\teffect_allele_frequency\tp_value\n"
    gz = _write_gz(tmp_path, header, [
        "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66",
    ], name="no_id.h.tsv.gz")
    out = tmp_path / "out.tsv"
    with pytest.raises(ValueError):
        harmonize_sumstats(gz, _study(), out)


def test_zero_rows_emits_warning(tmp_path, caplog):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tAC\tA\t2.0\tNA\t0.45\t0.66\tNA\tlo\t11\trs1\tNA\tNA",
    ], name="or_zero_rows.h.tsv.gz")
    out = tmp_path / "out.tsv"
    with caplog.at_level("WARNING"):
        n = harmonize_sumstats(gz, _study(), out)
    assert n == 0
    assert any(
        "GCST90091914" in record.message or "0 variants" in record.message
        for record in caplog.records
        if record.levelname == "WARNING"
    )


_BETA_CI_HEADER = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
                   "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t"
                   "hm_coordinate_conversion\thm_code\trsid\tci_upper\tci_lower\n")


def test_ci_derived_se_uses_the_linear_scale_for_beta_studies(tmp_path):
    """GWAS-SSF CI bounds sit on the effect's own scale. For a beta study that
    is the log-odds/linear scale already, so logging them is a scale error."""
    gz = _write_gz(tmp_path, _BETA_CI_HEADER, [
        "1\t10177\tA\tG\t0.10\tNA\t0.30\t1e-4\tNA\tlo\t11\trs1\t0.15\t0.05",
    ])
    out = tmp_path / "out.tsv"
    harmonize_sumstats(gz, _study(), out)
    got = pd.read_csv(out, sep="\t")["se"].iloc[0]
    assert got == pytest.approx((0.15 - 0.05) / (2 * _QNORM_975), rel=1e-6)


def test_ci_derived_se_stays_logarithmic_for_odds_ratio_studies(tmp_path):
    gz = _write_gz(tmp_path, _OR_HEADER, [
        "1\t10177\tA\tG\t1.10\tNA\t0.30\t1e-4\tNA\tlo\t11\trs1\t1.15\t1.05",
    ])
    out = tmp_path / "out.tsv"
    harmonize_sumstats(gz, _study(), out)
    got = pd.read_csv(out, sep="\t")["se"].iloc[0]
    expected = (math.log(1.15) - math.log(1.05)) / (2 * _QNORM_975)
    assert got == pytest.approx(expected, rel=1e-6)


def test_beta_scale_ci_spanning_zero_still_yields_an_se(tmp_path):
    """The `> 0` guard was written for OR bounds; a beta CI may be negative."""
    gz = _write_gz(tmp_path, _BETA_CI_HEADER, [
        "1\t10177\tA\tG\t-0.02\tNA\t0.30\t0.4\tNA\tlo\t11\trs1\t0.08\t-0.12",
    ])
    out = tmp_path / "out.tsv"
    harmonize_sumstats(gz, _study(), out)
    got = pd.read_csv(out, sep="\t")["se"].iloc[0]
    assert got == pytest.approx((0.08 - -0.12) / (2 * _QNORM_975), rel=1e-6)


def test_a_present_standard_error_is_never_overwritten_by_the_ci(tmp_path):
    gz = _write_gz(tmp_path, _BETA_CI_HEADER, [
        "1\t10177\tA\tG\t0.10\t0.0255\t0.30\t1e-4\tNA\tlo\t11\trs1\t0.15\t0.05",
    ])
    out = tmp_path / "out.tsv"
    harmonize_sumstats(gz, _study(), out)
    assert pd.read_csv(out, sep="\t")["se"].iloc[0] == pytest.approx(0.0255)


def test_resolve_rsid_coalesces_across_candidate_columns():
    """A partly-populated rsid column must not discard rows that hm_rsid covers."""
    df = pd.DataFrame({
        "rsid":       ["rs1", "", "NA"],
        "hm_rsid":    ["rsX", "rs2", ""],
        "variant_id": ["1_1_A_G", "1_2_A_G", "1_3_A_G"],
    })
    assert _resolve_rsid(df, "t.tsv").tolist() == ["rs1", "rs2", "1_3_A_G"]


def test_resolve_rsid_prefers_rsid_where_it_is_populated():
    df = pd.DataFrame({"rsid": ["rs1", "rs2"], "hm_rsid": ["rsX", "rsY"]})
    assert _resolve_rsid(df, "t.tsv").tolist() == ["rs1", "rs2"]


def test_resolve_rsid_returns_all_nan_when_every_candidate_is_empty():
    df = pd.DataFrame({"rsid": ["", "NA"], "hm_rsid": ["<NA>", "None"]})
    assert _resolve_rsid(df, "t.tsv").isna().all()


def test_resolve_rsid_raises_only_when_no_candidate_column_exists():
    with pytest.raises(ValueError, match="missing variant identifier"):
        _resolve_rsid(pd.DataFrame({"chromosome": [1]}), "t.tsv")


def test_sparse_rsid_column_no_longer_loses_the_rows_hm_rsid_covers(tmp_path):
    """End-to-end: 1 of 3 rows has a bare rsid, all 3 have hm_rsid."""
    header = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
              "standard_error\teffect_allele_frequency\tp_value\trsid\thm_rsid\n")
    gz = _write_gz(tmp_path, header, [
        "1\t100\tA\tG\t0.1\t0.02\t0.3\t1e-4\trs1\trs1",
        "1\t200\tA\tG\t0.2\t0.03\t0.3\t1e-4\t\trs2",
        "1\t300\tA\tG\t0.3\t0.04\t0.3\t1e-4\tNA\trs3",
    ])
    out = tmp_path / "out.tsv"
    n = harmonize_sumstats(gz, _study(), out)
    assert n == 3
    assert pd.read_csv(out, sep="\t")["rsid"].tolist() == ["rs1", "rs2", "rs3"]
