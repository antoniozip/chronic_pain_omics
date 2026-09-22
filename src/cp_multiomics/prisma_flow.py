"""Counts for the PRISMA 2020 flow diagram, read from the screening records.

The review's evidence comes from repositories, not from papers: every analysed
dataset was identified in GEO, the GWAS Catalog, PRIDE, MetaboLights or
Metabolomics Workbench and screened against its own deposited metadata and
files. The flow diagram therefore follows datasets through those repositories.
It used to follow the PubMed search instead, which supplied no analysed
dataset, and drew a PubMed/Europe PMC split that the plotting code invented by
halving the total (Europe PMC returned no records in the only search run).

Every count here is read from a tracked screening sheet, so the diagram cannot
state a number the records do not support. Exclusion reason codes are mapped to
display groups explicitly, and an unmapped code raises rather than falling into
an "other" bucket: a new reason must be given a place in the diagram on purpose.
Only once grouped are a source's smallest reasons folded together for drawing
(:meth:`SourceFlow.condensed`), and the fold is recorded in the sidecar.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import pandas as pd

PRISMA_DIR = Path("literature") / "prisma"
SUPP_DIR = Path("manuscript") / "supplementary"
META_DIR = Path("results") / "meta"

# Exclusion reasons drawn per repository before the rest fold into one line.
# Every reason still reaches the reader: Supplementary Tables 2, 6, 8 and 11
# carry each screened record with its reason code. The diagram cannot carry
# them all and print at a legible size: at 12 pt, 174 mm wide, the full
# breakdown made it 389 mm tall against a page of about 229 mm. Three reasons
# cost 20 mm more than two at any font size tried.
MAX_REASONS = 2
OTHER_REASONS = "Other reasons"


@dataclass(frozen=True)
class SourceFlow:
    """One repository's path from identification to inclusion."""

    source: str
    identified: int
    identified_label: str
    removed: tuple[tuple[str, int], ...] = field(default=())    # before screening
    excluded: tuple[tuple[str, int], ...] = field(default=())   # with reasons
    awaiting: tuple[tuple[str, int], ...] = field(default=())   # not yet decided
    included: int = 0
    included_label: str = "datasets"
    units: int | None = None                                    # study units pooled

    @property
    def n_excluded(self) -> int:
        return sum(n for _, n in self.excluded)

    @property
    def n_awaiting(self) -> int:
        return sum(n for _, n in self.awaiting)

    def condensed(self, keep: int = MAX_REASONS) -> SourceFlow:
        """The flow as drawn: its `keep` largest exclusion reasons, the rest folded."""
        return replace(self, excluded=top_reasons(self.excluded, keep))

    def check(self) -> None:
        """Identified must equal removed + excluded + awaiting + included."""
        total = (sum(n for _, n in self.removed) + self.n_excluded
                 + self.n_awaiting + self.included)
        if total != self.identified:
            raise ValueError(
                f"{self.source}: {self.identified} identified but "
                f"{total} accounted for; the flow does not add up")


# -- reason-code groups ------------------------------------------------------

NO_CONTRAST = "No pain-versus-control contrast"
NO_SPLIT = "No case-control split in the metadata"

GEO_EXCLUDED = {
    "no_case_control_split_in_series_matrix": NO_SPLIT,
    "no_pain_versus_control_contrast": NO_CONTRAST,
    "no_pain_phenotype": NO_CONTRAST,
    "in_vitro_hormone_treatment_not_patient_contrast": NO_CONTRAST,
    "in_vitro_stimulus_not_patient_contrast": NO_CONTRAST,
    "cultured_cells_across_passage_and_cell_type": NO_CONTRAST,
    "control_arm_is_pain_without_endometriosis": NO_CONTRAST,
    "superseries_duplicates_subseries": "SuperSeries duplicating pooled SubSeries",
    "single_cell_not_bulk_expression": "Single-cell deposit (own arm)",
    "herv_array_features_are_not_genes": "Features not genes, or deposit unreadable",
    "characteristics_fields_misaligned_in_deposit":
        "Features not genes, or deposit unreadable",
}
SINGLE_CELL_EXCLUDED = {
    "all_samples_case_like_no_control_arm": "No control arm or case-control split",
    "no_case_control_split_in_series_matrix": "No control arm or case-control split",
    "below_min_samples_per_group_1_case_3_control": "Too few samples per group",
    "bulk_rna_seq_not_single_cell_in_this_deposit": "Bulk, not single-cell, data",
}
GWAS_EXCLUDED = {
    "excluded: cohort/phenotype overlap (kept larger N)":
        "Cohort/phenotype overlap (larger kept)",
    "excluded: gene-based burden study": "Gene-based burden study",
    "excluded: no harmonized file available": "No harmonized summary statistics",
}
PRIDE_TRIAGE = {
    "ineligible": "Organism not eligible",
    "raw_only": "No processed quantification table",
}
PRIDE_REVIEW = {"no_pain_contrast": NO_CONTRAST}
MTBLS_REVIEW = {
    "no_pain_phenotype": "No pain phenotype",
    "no_pain_versus_control_contrast": NO_CONTRAST,
    "intervention_pre_post_not_case_control": NO_CONTRAST,
    "wrong_species": "Species not human, mouse or rat",
    "duplicate_of_MTBLS13513": "Duplicate deposit",
    "duplicate_of_MTBLS9662": "Duplicate deposit",
}
MTBLS_FILES = {
    "maf_columns_declared_but_empty": "Abundance columns empty",
    "no_control_arm_in_factors": "No control arm",
}
WORKBENCH_EXCLUDED = {
    "no_control_arm_in_factors": NO_CONTRAST,
    "no_pain_versus_control_contrast": NO_CONTRAST,
    "no_pain_phenotype": "No pain phenotype",
    "no_abundances_deposited": "No abundances deposited",
    "intervention_pre_post_not_case_control": "Other design or compartment",
    "faecal_scfa_panel_not_a_comparable_compartment": "Other design or compartment",
}


def top_reasons(items: tuple[tuple[str, int], ...],
                keep: int) -> tuple[tuple[str, int], ...]:
    """The `keep` largest reasons, largest first, the rest as one "other" line.

    A single leftover reason keeps its own name: folding it saves no line and
    would hide a reason for nothing.
    """
    ranked = sorted(items, key=lambda kv: (-kv[1], kv[0]))
    if len(ranked) <= keep + 1:
        return tuple(ranked)
    return tuple(ranked[:keep]) + ((OTHER_REASONS, sum(n for _, n in ranked[keep:])),)


def _group(codes: pd.Series, mapping: dict[str, str], source: str) -> tuple[tuple[str, int], ...]:
    """Count codes by display group, largest first; an unmapped code raises."""
    unknown = sorted(set(codes) - set(mapping))
    if unknown:
        raise ValueError(f"{source}: unmapped reason code(s) {unknown}; "
                         "add them to the mapping in prisma_flow.py")
    counts = codes.map(mapping).value_counts()
    return tuple((label, int(n)) for label, n in
                 sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _pooled_units(root: Path, modality: str) -> int | None:
    path = root / META_DIR / modality / "pooled_effects.csv"
    if not path.exists():
        return None
    ids = pd.read_csv(path, usecols=["study_ids"])["study_ids"].astype(str)
    return int(ids.str.split(";").explode().str.strip().nunique())


def geo_flow(root: Path) -> SourceFlow:
    d = pd.read_csv(root / PRISMA_DIR / "transcriptomics_candidates.csv")
    exc = d[d["verdict"] == "exclude"]
    pend = d[d["verdict"] == "pending"]
    # Pending is split by rule, not by table: a study cleared on design and
    # waiting only for its files is a different state from one held for a
    # design or file problem, whatever that problem's code.
    ready = int((pend["reason_code"] == "usable_pending_file_check").sum())
    units_path = root / SUPP_DIR / "Table_S1_transcriptomic_studies.csv"
    units = len(pd.read_csv(units_path)) if units_path.exists() else None
    return SourceFlow(
        source="GEO (bulk transcriptomics)", identified=len(d),
        identified_label="series",
        excluded=_group(exc["reason_code"], GEO_EXCLUDED, "GEO"),
        awaiting=(("Cleared on design; files not yet checked", ready),
                  ("Design or file problem under review", len(pend) - ready)),
        included=int((d["verdict"] == "include").sum()), included_label="accessions",
        units=units)


def single_cell_flow(root: Path) -> SourceFlow:
    d = pd.read_csv(root / PRISMA_DIR / "single_cell_candidates.csv")
    return SourceFlow(
        source="GEO single-cell series", identified=len(d),
        identified_label="series, reassessed",
        excluded=_group(d.loc[d["verdict"] == "exclude", "reason_code"],
                        SINGLE_CELL_EXCLUDED, "single-cell"),
        included=int((d["verdict"] == "include").sum()), included_label="series")


def gwas_flow(root: Path) -> SourceFlow:
    d = pd.read_csv(root / SUPP_DIR / "Table_S2_GWAS_studies.csv")
    return SourceFlow(
        source="GWAS Catalog", identified=len(d), identified_label="studies",
        excluded=_group(d.loc[d["Disposition"] == "Excluded", "Reason"],
                        GWAS_EXCLUDED, "GWAS"),
        included=int((d["Disposition"] == "Harmonized").sum()),
        included_label="studies")


def pride_flow(root: Path) -> SourceFlow:
    t = pd.read_csv(root / PRISMA_DIR / "proteomics_triage.csv")
    r = pd.read_csv(root / PRISMA_DIR / "proteomics_manual_review.csv")
    screened = t[t["verdict"] != "processed_table"]
    excluded = (_group(screened["verdict"], PRIDE_TRIAGE, "PRIDE")
                + _group(r.loc[r["verdict"] == "exclude", "reason_code"],
                         PRIDE_REVIEW, "PRIDE review"))
    if len(r) != int((t["verdict"] == "processed_table").sum()):
        raise ValueError("PRIDE: the manual review does not cover every processed table")
    return SourceFlow(
        source="PRIDE", identified=len(t), identified_label="datasets",
        excluded=excluded, included=int((r["verdict"] == "include").sum()),
        units=_pooled_units(root, "proteomics"))


def metabolights_flow(root: Path) -> SourceFlow:
    r = pd.read_csv(root / PRISMA_DIR / "metabolomics_manual_review.csv")
    f = pd.read_csv(root / PRISMA_DIR / "metabolomics_file_triage.csv")
    # Five screening exclusions carry reason_code "not_recorded" beside the
    # reason the reviewer recorded as proposed; that recorded reason is used.
    code = r["reason_code"].where(r["reason_code"] != "not_recorded",
                                  r["reason_code_proposed"])
    screened_out = _group(code[r["verdict"] == "exclude"], MTBLS_REVIEW, "MetaboLights")
    if len(f) != int((r["verdict"] == "include").sum()):
        raise ValueError("MetaboLights: file triage does not cover every study screened in")
    usable = f["verdict"] == "usable"
    return SourceFlow(
        source="MetaboLights", identified=len(r), identified_label="studies",
        excluded=screened_out + _group(f.loc[~usable, "reason_code"], MTBLS_FILES,
                                       "MetaboLights files"),
        included=int(usable.sum()), included_label="studies")


def workbench_flow(root: Path) -> SourceFlow:
    counts = pd.read_csv(root / PRISMA_DIR / "metabolomics_workbench_counts.csv")
    cands = pd.read_csv(root / PRISMA_DIR / "metabolomics_workbench_candidates.csv")
    enumerated = int(counts["n_enumerated"].iloc[0])
    removed = int(counts["excluded_no_pain_phenotype"].iloc[0])
    if enumerated - removed != len(cands):
        raise ValueError("Workbench: enumeration minus title screen is not the candidate count")
    return SourceFlow(
        source="Metabolomics Workbench", identified=enumerated,
        identified_label="studies enumerated",
        removed=(("Removed by automated title screen", removed),),
        excluded=_group(cands.loc[cands["verdict"] == "exclude", "reason_code"],
                        WORKBENCH_EXCLUDED, "Workbench"),
        included=int((cands["verdict"] == "include").sum()), included_label="studies")


def repository_flows(root: Path) -> list[SourceFlow]:
    """Every source's flow, each checked to add up."""
    flows = [geo_flow(root), single_cell_flow(root), gwas_flow(root),
             pride_flow(root), metabolights_flow(root), workbench_flow(root)]
    for flow in flows:
        flow.check()
    return flows


def metabolomic_units(root: Path) -> int | None:
    """Study units pooled across both metabolomic repositories."""
    return _pooled_units(root, "metabolomics")


def flow_table(flows: list[SourceFlow]) -> pd.DataFrame:
    """One row per drawn count: what the figure states, in long form."""
    rows: list[dict[str, object]] = []
    for f in flows:
        rows.append({"source": f.source, "stage": "identified", "label": "", "n": f.identified})
        for stage, items in (("removed", f.removed), ("excluded", f.excluded),
                             ("awaiting", f.awaiting)):
            rows += [{"source": f.source, "stage": stage, "label": label, "n": n}
                     for label, n in items]
        rows.append({"source": f.source, "stage": "included", "label": "", "n": f.included})
        if f.units is not None:
            rows.append({"source": f.source, "stage": "units", "label": "", "n": f.units})
    return pd.DataFrame(rows, columns=["source", "stage", "label", "n"])
