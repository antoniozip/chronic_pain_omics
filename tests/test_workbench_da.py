"""Contrast declaration and arm selection for the Workbench studies.

The recurring failure these guard against is a sample landing in the wrong
arm. A treated arm folded into the cases measures the treatment; folded into
the controls it measures protection. MTBLS5667's herbal-formula arm set that
trap, and ST000676 carries three arms with the same shape.
"""
from __future__ import annotations

import pytest

from cp_multiomics.metabolomics import workbench_da as wda
from cp_multiomics.metabolomics.workbench_mwtab import Analysis


def make(samples, factors):
    return Analysis(analysis_id="AN1", samples=samples, factors=factors,
                    values={"Taurine": [1.0] * len(samples)})


class TestFactorParsing:
    def test_multi_factor_cell(self):
        got = wda.parse_factor_cell(
            "AGE:16 week | DIET:HF STZ | SAMPLE_TYPE:Dorsal root ganglia")
        assert got == {"age": "16 week", "diet": "hf stz",
                       "sample_type": "dorsal root ganglia"}

    def test_case_and_whitespace_are_normalised(self):
        """ST000676 deposits both `Control` and `control` as diet levels; they
        are one arm and must merge."""
        assert (wda.parse_factor_cell("DIET:Control")["diet"]
                == wda.parse_factor_cell("diet:  control ")["diet"])

    def test_empty_and_malformed_cells(self):
        assert wda.parse_factor_cell("") == {}
        assert wda.parse_factor_cell("no colon here") == {}


class TestLevelMatching:
    def test_exact_levels(self):
        assert wda.level_matches("control", ["Control"])
        assert not wda.level_matches("hf stz", ["control"])

    def test_regex_spec(self):
        """ST003954 labels every case individually: IBS-C-1 ... IBS-C-61."""
        assert wda.level_matches("ibs-c-17", "re:^ibs")
        assert not wda.level_matches("hc", "re:^ibs")

    def test_bare_string_spec_is_refused(self):
        """A plain string would silently match nothing as a level list."""
        with pytest.raises(ValueError):
            wda.level_matches("control", "control")


class TestArmSelection:
    def test_declared_arms_only(self):
        """ST000676: HF, HF DR and HF STZ DR are named in neither arm and must
        be dropped, not assigned."""
        a = make(["s1", "s2", "s3", "s4", "s5"],
                 ["DIET:HF STZ", "DIET:Control", "DIET:HF",
                  "DIET:HF DR", "DIET:HF STZ DR"])
        spec = {"case": {"DIET": ["hf stz"]}, "control": {"DIET": ["control"]}}
        case, ctrl = wda.arm_columns(a, spec)
        assert case == ["s1"] and ctrl == ["s2"]

    def test_hf_stz_dr_does_not_match_hf_stz(self):
        """Exact level matching, not substring: the reversal arm is distinct."""
        assert not wda.level_matches("hf stz dr", ["hf stz"])

    def test_lean_arm_excluded_for_matched_control(self):
        """ST001412: comparing obese-neuropathy against lean would confound
        neuropathy with obesity."""
        a = make(["s1", "s2", "s3"],
                 ["Group:Obese neuropathy", "Group:Obese non neuropathy", "Group:Lean"])
        case, ctrl = wda.arm_columns(a, wda.WORKBENCH_STUDIES["ST001412"])
        assert case == ["s1"] and ctrl == ["s2"]

    def test_sample_missing_the_factor_joins_neither_arm(self):
        a = make(["s1", "s2"], ["Group:Control", "Stage:III"])
        spec = {"case": {"Group": ["endometriosis"]}, "control": {"Group": ["control"]}}
        case, ctrl = wda.arm_columns(a, spec)
        assert case == [] and ctrl == ["s1"]

    def test_all_declared_factors_must_match(self):
        a = make(["s1", "s2"], ["DIET:HF STZ | AGE:16 week", "DIET:HF STZ | AGE:24 week"])
        spec = {"case": {"DIET": ["hf stz"], "AGE": ["16 week"]},
                "control": {"DIET": ["control"]}}
        case, _ = wda.arm_columns(a, spec)
        assert case == ["s1"]

    def test_regex_arm_selection_end_to_end(self):
        a = make(["s1", "s2", "s3"], ["Group1:HC", "Group1:IBS-C-1", "Group1:IBS-D-7"])
        case, ctrl = wda.arm_columns(a, wda.WORKBENCH_STUDIES["ST003954"])
        assert case == ["s2", "s3"] and ctrl == ["s1"]


class TestSplitting:
    def test_split_restricts_both_arms_to_one_tissue(self):
        a = make(["p1", "p2", "d1", "d2"],
                 ["DIET:HF STZ | SAMPLE_TYPE:plasma",
                  "DIET:Control | SAMPLE_TYPE:plasma",
                  "DIET:HF STZ | SAMPLE_TYPE:Dorsal root ganglia",
                  "DIET:Control | SAMPLE_TYPE:Dorsal root ganglia"])
        spec = wda.WORKBENCH_STUDIES["ST000676"]
        case, ctrl = wda.arm_columns(a, spec, split_value="plasma")
        assert case == ["p1"] and ctrl == ["p2"]
        case, ctrl = wda.arm_columns(a, spec, split_value="dorsal root ganglia")
        assert case == ["d1"] and ctrl == ["d2"]

    def test_split_values_are_discovered_in_order(self):
        a = make(["a", "b", "c"], ["SAMPLE_TYPE:plasma", "SAMPLE_TYPE:Sciatic nerve",
                                   "SAMPLE_TYPE:plasma"])
        assert wda.split_values(a, "SAMPLE_TYPE") == ["plasma", "sciatic nerve"]

    def test_derived_ids_follow_the_gse241361_shape(self):
        """Species and parentage are recovered by suffix-stripping, so the
        separator must not appear inside the accession."""
        assert wda.derived_study_id("ST000676", "dorsal root ganglia") == \
            "ST000676_Dorsal_root_ganglia"
        assert wda.derived_study_id("ST003984", "peritoneal fluid") == \
            "ST003984_Peritoneal_fluid"
        assert wda.derived_study_id("ST000676", "plasma").rsplit("_", 1)[0] == "ST000676"


class TestFrame:
    def test_refmet_name_preferred_over_vendor_name(self):
        a = Analysis(analysis_id="AN1", samples=["s1", "s2"],
                     factors=["Group:Control", "Group:Case"],
                     values={"Ac-Orn": [1.0, 2.0], "ADMA": [3.0, 4.0]},
                     refmet={"Ac-Orn": "N-Acetylornithine"})
        frame, names = wda.analysis_frame(a)
        assert list(names) == ["N-Acetylornithine", "ADMA"]
        assert list(frame.columns) == ["s1", "s2"]
        assert frame.iloc[0].tolist() == [1.0, 2.0]


class TestSpecsAreWellFormed:
    @pytest.mark.parametrize("acc", sorted(wda.WORKBENCH_STUDIES))
    def test_every_spec_declares_both_arms_and_a_condition(self, acc):
        spec = wda.WORKBENCH_STUDIES[acc]
        assert spec["case"] and spec["control"]
        assert spec["condition"]
        # A split study takes its region from the split value; every other
        # study must name one, since `region` is a column in the effects table.
        assert ("region" in spec) != ("split_by" in spec)

    @pytest.mark.parametrize("acc", sorted(wda.WORKBENCH_STUDIES))
    def test_arms_are_disjoint(self, acc):
        """A level in both arms would put the same samples on both sides."""
        spec = wda.WORKBENCH_STUDIES[acc]
        for key, case_spec in spec["case"].items():
            ctrl_spec = spec["control"].get(key)
            if isinstance(case_spec, list) and isinstance(ctrl_spec, list):
                assert not ({wda.norm_level(s) for s in case_spec}
                            & {wda.norm_level(s) for s in ctrl_spec})
            elif isinstance(case_spec, str) and isinstance(ctrl_spec, list):
                for level in ctrl_spec:
                    assert not wda.level_matches(wda.norm_level(level), case_spec)
