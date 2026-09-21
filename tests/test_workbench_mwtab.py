"""Parsing Metabolomics Workbench mwtab files.

The tests that matter most here are about what a cell means when it is not a
number. ST003984 deposits 6,300 `&lt; LOD` cells and ST003177 41,245 `NA`;
reading any of those as zero would manufacture effect sizes in precisely the
case-versus-control comparisons the arm exists to make.
"""
from __future__ import annotations

import math

import pytest

from cp_multiomics.metabolomics import workbench_mwtab as mw

# Two analyses over the same samples — the ion-mode pair that must never be
# treated as two studies.
TWO_ANALYSES = "\n".join([
    "#METABOLOMICS WORKBENCH someone_1 DATATRACK_ID:1 STUDY_ID:ST000001 ANALYSIS_ID:AN000001",
    "#MS_METABOLITE_DATA",
    "MS_METABOLITE_DATA:UNITS\tng/ml",
    "MS_METABOLITE_DATA_START",
    "Samples\tS1\tS2\tS3\tS4",
    "Factors\tGroup:Control\tGroup:Control\tGroup:Case\tGroup:Case",
    "Taurine\t1.0\t2.0\t3.0\t4.0",
    "Creatine\t10\t20\t30\t40",
    "MS_METABOLITE_DATA_END",
    "#METABOLITES",
    "METABOLITES_START",
    "metabolite_name\tRefMet Name\tWorkBench Metabolite_ID",
    "Taurine\tTaurine\tME1",
    "Creatine\tCreatine\tME2",
    "METABOLITES_END",
    "#END",
    "#METABOLOMICS WORKBENCH someone_1 DATATRACK_ID:1 STUDY_ID:ST000001 ANALYSIS_ID:AN000002",
    "MS_METABOLITE_DATA_START",
    "Samples\tS1\tS2\tS3\tS4",
    "Factors\tGroup:Control\tGroup:Control\tGroup:Case\tGroup:Case",
    "Betaine\t5\t6\t7\t8",
    "MS_METABOLITE_DATA_END",
    "#END",
])


class TestParseValue:
    @pytest.mark.parametrize("raw", [
        "&lt; LOD",     # HTML-escaped, as deposited
        "&lt; LLOQ",
        "< LOD",
        "<LLOQ",
        "NA", "na", "N/A", "ND", "NaN", "null",
        "", "   ", "-", ".",
        "below detection limit",
    ])
    def test_non_measurements_are_nan_never_zero(self, raw):
        """A below-detection marker means "not measurable", not "abundance 0".

        Substituting zero for a compound absent from every control and present
        in every case yields an enormous fabricated effect size.
        """
        value = mw.parse_value(raw)
        assert math.isnan(value), raw
        assert value != 0.0

    @pytest.mark.parametrize("raw,expected", [
        ("127.443", 127.443), ("0", 0.0), ("-1.5", -1.5),
        ("1e3", 1000.0), (" 42 ", 42.0), ("1,234.5", 1234.5),
    ])
    def test_real_numbers_parse(self, raw, expected):
        assert mw.parse_value(raw) == pytest.approx(expected)

    def test_a_genuine_zero_is_kept(self):
        """Only markers become NaN; a deposited 0 is a measurement."""
        assert mw.parse_value("0.0") == 0.0
        assert not math.isnan(mw.parse_value("0.0"))

    def test_none_is_nan(self):
        assert math.isnan(mw.parse_value(None))


class TestSplitAnalyses:
    def test_both_analyses_are_found(self):
        pairs = mw.split_analyses(TWO_ANALYSES)
        assert [a for a, _ in pairs] == ["AN000001", "AN000002"]

    def test_each_document_holds_only_its_own_data(self):
        pairs = dict(mw.split_analyses(TWO_ANALYSES))
        assert "Taurine" in pairs["AN000001"]
        assert "Taurine" not in pairs["AN000002"]
        assert "Betaine" in pairs["AN000002"]

    def test_text_without_a_header_yields_nothing(self):
        """A truncated download must not be read as one untitled analysis."""
        assert mw.split_analyses("Samples\tS1\nTaurine\t1.0") == []
        assert mw.split_analyses("") == []


class TestParseAnalysis:
    def test_samples_factors_and_values(self):
        a = mw.parse_mwtab(TWO_ANALYSES)[0]
        assert a.analysis_id == "AN000001"
        assert a.samples == ["S1", "S2", "S3", "S4"]
        assert a.factors[0] == "Group:Control" and a.factors[3] == "Group:Case"
        assert a.values["Taurine"] == [1.0, 2.0, 3.0, 4.0]
        assert a.units == "ng/ml"
        assert a.n_samples == 4

    def test_two_analyses_are_separate_units(self):
        """Same samples, two ion modes. Collapsing happens downstream; here
        they must stay distinct so neither is lost."""
        analyses = mw.parse_mwtab(TWO_ANALYSES)
        assert len(analyses) == 2
        assert set(analyses[0].values) == {"Taurine", "Creatine"}
        assert set(analyses[1].values) == {"Betaine"}

    def test_refmet_names_are_read(self):
        a = mw.parse_mwtab(TWO_ANALYSES)[0]
        assert a.refmet["Taurine"] == "Taurine"

    def test_refmet_absent_is_empty_not_an_error(self):
        a = mw.parse_mwtab(TWO_ANALYSES)[1]
        assert a.refmet == {}

    def test_refmet_dash_is_not_a_name(self):
        doc = TWO_ANALYSES.replace("Taurine\tTaurine\tME1", "Taurine\t-\tME1")
        a = mw.parse_mwtab(doc)[0]
        assert "Taurine" not in a.refmet

    def test_below_detection_row_keeps_its_measured_values(self):
        doc = TWO_ANALYSES.replace("Taurine\t1.0\t2.0\t3.0\t4.0",
                                   "Taurine\t&lt; LOD\t2.0\t&lt; LLOQ\t4.0")
        a = mw.parse_mwtab(doc)[0]
        row = a.values["Taurine"]
        assert math.isnan(row[0]) and row[1] == 2.0
        assert math.isnan(row[2]) and row[3] == 4.0

    def test_an_all_missing_metabolite_is_dropped(self):
        """A row with no measurement anywhere cannot yield an effect size."""
        doc = TWO_ANALYSES.replace("Creatine\t10\t20\t30\t40",
                                   "Creatine\t&lt; LOD\tNA\t&lt; LOD\t")
        a = mw.parse_mwtab(doc)[0]
        assert "Creatine" not in a.values
        assert "Taurine" in a.values

    def test_short_factor_row_is_padded_not_zipped(self):
        """Zipping a short Factors row against the samples would silently
        assign every later column to the wrong arm."""
        doc = TWO_ANALYSES.replace(
            "Factors\tGroup:Control\tGroup:Control\tGroup:Case\tGroup:Case",
            "Factors\tGroup:Control\tGroup:Control")
        a = mw.parse_mwtab(doc)[0]
        assert len(a.factors) == len(a.samples) == 4
        assert a.factors[2] == "" and a.factors[3] == ""

    def test_short_value_row_is_padded(self):
        doc = TWO_ANALYSES.replace("Taurine\t1.0\t2.0\t3.0\t4.0", "Taurine\t1.0\t2.0")
        a = mw.parse_mwtab(doc)[0]
        assert len(a.values["Taurine"]) == 4
        assert math.isnan(a.values["Taurine"][3])

    def test_missing_data_block_yields_no_analysis(self):
        doc = "\n".join([
            "#METABOLOMICS WORKBENCH x DATATRACK_ID:1 STUDY_ID:ST1 ANALYSIS_ID:AN000009",
            "#STUDY", "ST:STUDY_TITLE\tNo data here", "#END",
        ])
        assert mw.parse_mwtab(doc) == []

    def test_block_without_a_samples_row_is_refused(self):
        doc = TWO_ANALYSES.replace("Samples\tS1\tS2\tS3\tS4", "Bins\tS1\tS2\tS3\tS4", 1)
        analyses = mw.parse_mwtab(doc)
        assert [a.analysis_id for a in analyses] == ["AN000002"]
