"""Tests for the single-cell arm's group assignment and pseudobulk readers.

Each of these covers a defect that produced output of exactly the right shape
and the wrong content, which is the only kind this arm can suffer from: there
is no per-sample ground truth to compare against, so a matrix summed along the
wrong axis, aligned on the wrong index, or missing a repeated animal reads as
a perfectly ordinary result.

  * **Orientation.** Three of the deposits write cells by genes and the rest
    genes by cells, and only a filename or a comment says which. Summing the
    wrong axis yields a vector of the right length whenever the matrix is
    close to square.
  * **Feature alignment.** GSE198608's two submission batches deposit the same
    32,883 genes in two different orders and two different spellings, so a
    positional merge transposes the whole annotation and an unmapped merge
    silently drops 2,816 genes as unmeasured.
  * **Duplicate symbols.** A 10x features file routinely gives two Ensembl ids
    one symbol. Keeping one occurrence discards the other's counts.
  * **Repeated animals.** GSE162807 sequenced three mice twice. Entering both
    libraries counts the animal twice, which is the pseudoreplication the
    whole arm exists to avoid at the cell level and would reintroduce at the
    animal level.
  * **Cell calling.** GSE254360 deposits the raw droplet matrix, 2.4 million
    barcodes for a few thousand cells. Summing every droplet adds the ambient
    profile to the sample.
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from cp_multiomics.single_cell import (  # noqa: E402
    Pseudobulk,
    read_dense_split,
    read_lines,
    read_mtx,
    strip_version,
)


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, REPO_ROOT / "pipeline" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


step18 = _load("pipeline_18", "18_assign_single_cell_groups.py")
step19 = _load("pipeline_19", "19_pseudobulk_single_cell.py")


# ---------------------------------------------------------------------------
# Matrix Market orientation
# ---------------------------------------------------------------------------

def _write_mtx(path: Path, rows: int, cols: int, entries: list[tuple[int, int, int]]):
    with gzip.open(path, "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("%a comment\n")
        handle.write(f"{rows} {cols} {len(entries)}\n")
        for r, c, v in entries:
            handle.write(f"{r} {c} {v}\n")


def test_mtx_genes_by_cells_sums_over_cells(tmp_path: Path) -> None:
    # 3 genes x 2 cells. Gene 1 has 5+7, gene 2 has 0, gene 3 has 11.
    path = tmp_path / "m.mtx.gz"
    _write_mtx(path, 3, 2, [(1, 1, 5), (1, 2, 7), (3, 2, 11)])
    block = read_mtx(path, ["A", "B", "C"], "S1", min_umi=0, cells_are_rows=False)
    assert list(block.columns["S1"]) == [12.0, 0.0, 11.0]
    assert block.cells["S1"] == (2, 2)


def test_mtx_cells_by_genes_sums_over_the_other_axis(tmp_path: Path) -> None:
    # The same data transposed: 2 cells x 3 genes must give the same totals.
    path = tmp_path / "m.mtx.gz"
    _write_mtx(path, 2, 3, [(1, 1, 5), (2, 1, 7), (2, 3, 11)])
    block = read_mtx(path, ["A", "B", "C"], "S1", min_umi=0, cells_are_rows=True)
    assert list(block.columns["S1"]) == [12.0, 0.0, 11.0]
    assert block.cells["S1"] == (2, 2)


def test_mtx_refuses_a_feature_list_of_the_wrong_length(tmp_path: Path) -> None:
    path = tmp_path / "m.mtx.gz"
    _write_mtx(path, 3, 2, [(1, 1, 5)])
    with pytest.raises(ValueError, match="gene rows against"):
        read_mtx(path, ["A", "B"], "S1", min_umi=0, cells_are_rows=False)


def test_umi_floor_drops_a_barcode_below_it(tmp_path: Path) -> None:
    # Cell 2 carries 3 UMIs and is an empty droplet at a floor of 10.
    path = tmp_path / "m.mtx.gz"
    _write_mtx(path, 2, 2, [(1, 1, 40), (2, 1, 10), (1, 2, 3)])
    kept = read_mtx(path, ["A", "B"], "S1", min_umi=10, cells_are_rows=False)
    assert list(kept.columns["S1"]) == [40.0, 10.0]
    assert kept.cells["S1"] == (2, 1)
    everything = read_mtx(path, ["A", "B"], "S1", min_umi=0, cells_are_rows=False)
    assert list(everything.columns["S1"]) == [43.0, 10.0]


# ---------------------------------------------------------------------------
# Dense tables, split by a prefix on the cell name
# ---------------------------------------------------------------------------

def _write_dense(path: Path, columns: list[str], rows: dict[str, list[int]]):
    frame = pd.DataFrame(rows, index=columns).T
    frame.index.name = ""
    frame.to_csv(path, sep="\t", compression="gzip")


def test_dense_split_sums_within_each_sample(tmp_path: Path) -> None:
    path = tmp_path / "t.txt.gz"
    _write_dense(path, ["SI01_AAA", "SI01_BBB", "SNI01_CCC"],
                 {"Gene1": [1, 2, 30], "Gene2": [4, 5, 60]})
    block = read_dense_split(
        path, lambda c: c.split("_")[0] if c.startswith(("SI", "SNI")) else None,
        min_umi=0)
    assert block.features == ["Gene1", "Gene2"]
    assert list(block.columns["SI01"]) == [3.0, 9.0]
    assert list(block.columns["SNI01"]) == [30.0, 60.0]
    assert block.cells["SI01"] == (2, 2)


def test_dense_split_refuses_a_table_no_column_maps_into(tmp_path: Path) -> None:
    path = tmp_path / "t.txt.gz"
    _write_dense(path, ["X_1"], {"Gene1": [1]})
    with pytest.raises(ValueError, match="no column mapped"):
        read_dense_split(path, lambda _c: None, min_umi=0)


# ---------------------------------------------------------------------------
# Feature alignment across samples
# ---------------------------------------------------------------------------

def test_align_features_matches_on_id_not_position() -> None:
    a = Pseudobulk(["A", "B", "C"], {"S1": np.array([1.0, 2.0, 3.0])}, {"S1": (1, 1)})
    b = Pseudobulk(["C", "B", "A"], {"S2": np.array([30.0, 20.0, 10.0])}, {"S2": (1, 1)})
    merged = step19.align_features("GSEX", [a, b])
    assert merged.features == ["A", "B", "C"]
    assert list(merged.columns["S1"]) == [1.0, 2.0, 3.0]
    assert list(merged.columns["S2"]) == [10.0, 20.0, 30.0]


def test_align_features_drops_what_is_not_measured_everywhere() -> None:
    a = Pseudobulk(["A", "B"], {"S1": np.array([1.0, 2.0])}, {"S1": (1, 1)})
    b = Pseudobulk(["B", "C"], {"S2": np.array([20.0, 30.0])}, {"S2": (1, 1)})
    merged = step19.align_features("GSEX", [a, b])
    # B only: filling A and C with zero would fabricate the largest effect in
    # the study out of an annotation difference.
    assert merged.features == ["B"]
    assert list(merged.columns["S1"]) == [2.0]
    assert list(merged.columns["S2"]) == [20.0]


def test_align_features_sums_duplicate_symbols_before_aligning() -> None:
    # Two Ensembl ids, one symbol -- the ordinary 10x features file.
    a = Pseudobulk(["A", "A", "B"], {"S1": np.array([1.0, 2.0, 5.0])}, {"S1": (1, 1)})
    b = Pseudobulk(["A", "B"], {"S2": np.array([10.0, 50.0])}, {"S2": (1, 1)})
    merged = step19.align_features("GSEX", [a, b])
    assert merged.features == ["A", "B"]
    assert list(merged.columns["S1"]) == [3.0, 5.0]


def test_align_features_refuses_a_study_with_no_shared_feature() -> None:
    a = Pseudobulk(["A"], {"S1": np.array([1.0])}, {"S1": (1, 1)})
    b = Pseudobulk(["B"], {"S2": np.array([1.0])}, {"S2": (1, 1)})
    with pytest.raises(ValueError, match="no feature is present in every sample"):
        step19.align_features("GSEX", [a, b])


# ---------------------------------------------------------------------------
# Repeated animals
# ---------------------------------------------------------------------------

def test_collapse_subjects_sums_two_libraries_of_one_animal() -> None:
    block = Pseudobulk(["A", "B"],
                       {"GSM1": np.array([1.0, 2.0]), "GSM2": np.array([10.0, 20.0]),
                        "GSM3": np.array([100.0, 200.0])},
                       {"GSM1": (5, 5), "GSM2": (7, 7), "GSM3": (9, 9)})
    groups = pd.DataFrame({"sample_id": ["GSM1", "GSM2", "GSM3"],
                           "subject": ["GSE_103", "GSE_103", "GSE_104"],
                           "group": ["case", "case", "control"]})
    merged = step19.collapse_subjects(block, groups)
    assert set(merged.columns) == {"GSE_103", "GSE_104"}
    assert list(merged.columns["GSE_103"]) == [11.0, 22.0]
    assert merged.cells["GSE_103"] == (12, 12)


def test_collapse_subjects_leaves_unrepeated_samples_alone() -> None:
    block = Pseudobulk(["A"], {"GSM1": np.array([1.0])}, {"GSM1": (5, 5)})
    groups = pd.DataFrame({"sample_id": ["GSM1"], "subject": ["GSM1"],
                           "group": ["case"]})
    merged = step19.collapse_subjects(block, groups)
    assert list(merged.columns) == ["GSM1"]


# ---------------------------------------------------------------------------
# Sidecar readers
# ---------------------------------------------------------------------------

def test_read_lines_takes_the_named_column(tmp_path: Path) -> None:
    path = tmp_path / "features.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("ENSG1\tAAA\tGene Expression\nENSG2\tBBB\tGene Expression\n")
    assert read_lines(path, column=1) == ["AAA", "BBB"]
    assert read_lines(path, column=0) == ["ENSG1", "ENSG2"]


def test_read_lines_skips_a_header_when_told_to(tmp_path: Path) -> None:
    path = tmp_path / "all_genes.csv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("gene_id,gene_name,genome\nENSRNOG1,Arsj,rat\n")
    assert read_lines(path, column=1, sep=",", skip_header=True) == ["Arsj"]


def test_strip_version_removes_the_ensembl_suffix() -> None:
    assert strip_version(["ENSMUSG00000102693.1", "ENSMUSG00000064842"]) == \
        ["ENSMUSG00000102693", "ENSMUSG00000064842"]


# ---------------------------------------------------------------------------
# Group assignment
# ---------------------------------------------------------------------------

def test_every_registry_study_is_either_ruled_or_recorded_as_excluded() -> None:
    registry = pd.read_csv(REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv")
    covered = set(step18.CONTRASTS) | set(step18.NO_CONTRAST)
    assert set(registry["accession"]) == covered, (
        "a study in the registry with neither a contrast rule nor an exclusion "
        "reason is dropped from the arm without a record of why")


def test_no_study_is_both_ruled_and_excluded() -> None:
    assert not (set(step18.CONTRASTS) & set(step18.NO_CONTRAST))


def test_every_contrast_rule_carries_a_note() -> None:
    for accession, rule in step18.CONTRASTS.items():
        assert rule.get("note"), f"{accession} has no note explaining its rule"
        assert "case" in rule and "control" in rule


def test_a_treatment_arm_is_not_matched_as_a_case() -> None:
    # GSE328175's SNI+Transplant animals are treated, and entering them as
    # cases is what the 2026-08-29 audit disqualified GSE343056 for.
    rule = step18.CONTRASTS["GSE328175"]
    treated = {"treatment": "SNI+Transplant", "title": "SNI_Transplant1"}
    assert not step18.matches(treated, rule["case"])
    assert not step18.matches(treated, rule["control"])
    assert step18.matches({"treatment": "SNI"}, rule["case"])


def test_superficial_injury_is_read_as_the_control_arm() -> None:
    # The keyword audit reported zero control-like fields for GSE134003
    # because "SI" is in no control vocabulary.
    rule = step18.CONTRASTS["GSE134003"]
    assert step18.matches({"title": "SI07"}, rule["control"])
    assert step18.matches({"title": "SNI07"}, rule["case"])
    # And the two must not both match: SNI starts with S, and a looser control
    # pattern would put every case in the control arm as well.
    assert not step18.matches({"title": "SNI07"}, rule["control"])


def test_a_study_below_the_minimum_per_arm_is_excluded() -> None:
    rows = [{"group": "case", "subject": "a"}, {"group": "control", "subject": "b"},
            {"group": "control", "subject": "c"}]
    verdict, reason, n_case, n_control = step18.study_verdict(rows)
    assert verdict == "exclude"
    assert (n_case, n_control) == (1, 2)
    assert "below_min_samples_per_group" in reason


def test_repeated_subjects_count_once_towards_the_arm_minimum() -> None:
    # Two libraries of one animal are one sample, so this is 1 v 2, not 2 v 2.
    rows = [{"group": "case", "subject": "x"}, {"group": "case", "subject": "x"},
            {"group": "control", "subject": "y"}, {"group": "control", "subject": "z"}]
    verdict, _reason, n_case, n_control = step18.study_verdict(rows)
    assert (n_case, n_control) == (1, 2)
    assert verdict == "exclude"


# ---------------------------------------------------------------------------
# The arm must not reach the bulk pool
# ---------------------------------------------------------------------------

def test_single_cell_is_not_one_of_step_05s_all_modalities() -> None:
    source = (REPO_ROOT / "pipeline" / "05_per_study_da.R").read_text()
    line = next(ln for ln in source.splitlines()
                if ln.startswith("ALL_MODALITIES <-"))
    assert "single_cell" not in line, (
        "putting single_cell in ALL_MODALITIES lets `--modality all` merge the "
        "pseudobulk arm into the bulk transcriptomic pool")


def test_single_cell_is_not_one_of_step_06s_all_modalities() -> None:
    source = (REPO_ROOT / "pipeline" / "06_meta_analysis.R").read_text()
    line = next(ln for ln in source.splitlines()
                if ln.startswith("ALL_MODALITIES <-"))
    assert "single_cell" not in line


def test_the_arm_writes_no_groups_file_into_the_shared_geo_cache() -> None:
    # Every one of these accessions is also in the bulk arm's screening sheet,
    # so a groups file left in data/raw/geo_cache/ is a live wire: it is what
    # the bulk transcriptomics DA reads.
    source = (REPO_ROOT / "pipeline" / "18_assign_single_cell_groups.py").read_text()
    assert 'GROUPS_DIR = REPO_ROOT / "data" / "interim" / "single_cell" / "groups"' \
        in source
    assert "geo_cache" in source, "the cache is still read for series matrices"
    assert 'CACHE / f"{accession}_groups.csv"' not in source


def test_no_study_is_included_in_both_the_single_cell_and_the_bulk_arm() -> None:
    # All ten of these accessions are also rows in the bulk arm's screening
    # sheet, excluded. Flipping one to `include` there would put the same
    # samples into two pools under two different effect tables -- the
    # double-counting `conf/analysis/superseded_studies.csv` exists to prevent,
    # arriving by a route that registry does not cover.
    prisma = REPO_ROOT / "literature" / "prisma"
    arm = pd.read_csv(prisma / "single_cell_candidates.csv")
    bulk = pd.read_csv(prisma / "transcriptomics_candidates.csv")
    in_arm = set(arm.loc[arm["verdict"] == "include", "accession"])
    in_bulk = set(bulk.loc[bulk["verdict"] == "include", "accession"])
    both = sorted(in_arm & in_bulk)
    assert not both, (
        f"{both} would enter the pseudobulk arm and the bulk transcriptomic "
        "pool on the same samples")


def test_the_arm_is_the_ten_studies_the_adjudication_settled_on() -> None:
    arm = pd.read_csv(
        REPO_ROOT / "literature" / "prisma" / "single_cell_candidates.csv")
    assert len(arm) == 18, "the screening record must cover all 18 retrieved series"
    assert (arm["verdict"] == "include").sum() == 10
    # Every exclusion carries a reason; a blank one is a study dropped silently.
    assert arm.loc[arm["verdict"] == "exclude", "reason_code"].notna().all()


def test_manifest_records_cells_for_a_repeated_animals_gsms() -> None:
    # The counts frame is keyed by subject and the cell counts by GSM. Reading
    # either through the other's key returns nothing and writes a zero, which
    # looks exactly like a sample that contributed no cells. GSE162807 recorded
    # 0 of 0 barcodes across all 28 of its samples that way.
    groups = pd.DataFrame({
        "sample_id": ["GSM1", "GSM2", "GSM3"],
        "subject": ["GSE_103", "GSE_103", "GSE_104"],
        "group": ["case", "case", "control"]})
    frame = pd.DataFrame({"GSE_103": [11, 22], "GSE_104": [100, 200]},
                         index=["A", "B"])
    cells = {"GSM1": (5, 4), "GSM2": (7, 6), "GSM3": (9, 9)}
    rows = step19.manifest_rows("GSEX", groups, cells, frame, missing=[])
    by_gsm = {r["sample_id"]: r for r in rows}
    assert by_gsm["GSM1"]["cells_seen"] == 5 and by_gsm["GSM1"]["cells_kept"] == 4
    assert by_gsm["GSM2"]["cells_seen"] == 7 and by_gsm["GSM2"]["cells_kept"] == 6
    # Both libraries of one animal report that animal's summed counts.
    assert by_gsm["GSM1"]["total_counts"] == by_gsm["GSM2"]["total_counts"] == 33
    assert by_gsm["GSM3"]["total_counts"] == 300


def test_manifest_marks_an_assigned_sample_absent_from_the_deposit() -> None:
    groups = pd.DataFrame({"sample_id": ["GSM1", "GSM2"], "subject": ["GSM1", "GSM2"],
                           "group": ["case", "control"]})
    frame = pd.DataFrame({"GSM1": [3]}, index=["A"])
    rows = step19.manifest_rows("GSEX", groups, {"GSM1": (2, 2)}, frame,
                                missing=["GSM2"])
    by_gsm = {r["sample_id"]: r for r in rows}
    assert by_gsm["GSM2"]["status"] == "absent_from_deposit"
    assert by_gsm["GSM2"]["group"] == "control"   # recorded, never dropped
    assert by_gsm["GSM1"]["status"] == "pseudobulked"


def test_a_merged_table_column_that_maps_to_nothing_is_refused() -> None:
    # A merged deposit holds that study's cells and nothing else, so a column
    # naming a library the map does not know is a broken map, not a missing
    # sample -- and the two look identical downstream. GSE155622's counts
    # tables spell a library `SNI.6h` where its metadata spells it `SNI 6h`,
    # so seven of ten libraries mapped to nothing and the study came out with
    # three controls and an empty case arm.
    from collections import Counter
    with pytest.raises(ValueError, match="the sample map does not know"):
        step19.assert_every_prefix_maps("GSEX", Counter({"SNI.6h": 7520}), ())


def test_a_prefix_declared_outside_the_arm_is_allowed() -> None:
    from collections import Counter
    step19.assert_every_prefix_maps("GSEX", Counter({"Excluded_1": 10}),
                                    ("Excluded_1",))


def test_the_gse155622_map_carries_both_spellings_the_deposit_uses() -> None:
    for name in ("SNI 6h", "SNI 24h_2", "SNI 14d"):
        dotted = name.replace(" ", ".")
        assert step19.GSE155622_PREFIXES[name] == step19.GSE155622_PREFIXES[dotted]
    # Ten libraries, seven of which have a space to double.
    assert len(set(step19.GSE155622_PREFIXES.values())) == 10
