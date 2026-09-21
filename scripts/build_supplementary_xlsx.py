#!/usr/bin/env python
"""Build the numbered supplementary workbooks the journal receives.

The manuscript used to cite its underlying data by repository path --
`conf/analysis/superseded_studies.csv` and a dozen others. A reader with the
PDF cannot open any of them, and will not clone a repository to see why a
study was excluded. Every such file is now a numbered Supplementary Table,
shipped as .xlsx beside the paper, and the text cites the number.

The CSVs stay where they are and stay authoritative: the pipeline writes them,
`tests/test_manuscript_claims.py` reads several of them, and this script only
renders them. Treat the workbooks as build artifacts, like the .docx files --
regenerate rather than edit, or the next build silently reverts the edit.

A table drawing on more than one CSV gets one sheet per source, named for the
file it came from, because merging them would lose which record said what.

Usage:
    python scripts/build_supplementary_xlsx.py
    python scripts/build_supplementary_xlsx.py --only 12 --only 14
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "manuscript" / "supplementary"

FONT = "Arial"
MAX_COL_WIDTH = 60          # a free-text column would otherwise run off the page
HEADER_FILL_ROWS = 1


@dataclass(frozen=True)
class Table:
    number: int
    slug: str                       # file-name stem after Table_SN_
    title: str
    sources: tuple[Path, ...]
    sheets: tuple[str, ...] = field(default=())

    def sheet_names(self) -> tuple[str, ...]:
        if self.sheets:
            return self.sheets
        return (self.slug[:31],)


def _p(rel: str) -> Path:
    return REPO_ROOT / rel


#: Every numbered supplementary data table, in order. S1-S3 already shipped as
#: CSV and are rendered here too, so the bundle is one format throughout.
TABLES: tuple[Table, ...] = (
    Table(1, "transcriptomic_studies", "Transcriptomic study units",
          (_p("manuscript/supplementary/Table_S1_transcriptomic_studies.csv"),)),
    Table(2, "GWAS_studies", "Retrieved GWAS studies",
          (_p("manuscript/supplementary/Table_S2_GWAS_studies.csv"),)),
    Table(3, "MAGMA_gene_results", "MAGMA gene-based results",
          (_p("manuscript/supplementary/Table_S3_MAGMA_gene_results.csv"),)),
    Table(4, "egger_test", "Egger's regression test per feature",
          (_p("results/meta/transcriptomics/egger_test.csv"),)),
    Table(5, "PXD013362_peptides", "Per-peptide proteomic effect sizes",
          (_p("results/meta/proteomics/PXD013362_pooled.csv"),)),
    Table(6, "proteomic_screening", "Proteomic screening and triage",
          (_p("literature/prisma/proteomics_triage.csv"),
           _p("literature/prisma/proteomics_manual_review.csv")),
          ("triage", "manual_review")),
    Table(7, "proteomic_pooled", "Pooled proteomic effect sizes",
          (_p("results/meta/proteomics/pooled_effects.csv"),)),
    Table(8, "metabolomic_screening", "Metabolomic screening and triage",
          (_p("literature/prisma/metabolomics_manual_review.csv"),
           _p("literature/prisma/metabolomics_file_triage.csv"),
           _p("literature/prisma/metabolomics_workbench_candidates.csv"),
           _p("literature/prisma/metabolomics_workbench_counts.csv"),
           _p("literature/prisma/metabolomics_workbench_validation.csv")),
          ("metabolights_review", "metabolights_files", "workbench_candidates",
           "workbench_counts", "workbench_validation")),
    Table(9, "metabolomic_pooled", "Pooled metabolomic effect sizes",
          (_p("results/meta/metabolomics/pooled_effects.csv"),)),
    Table(10, "contrast_audit", "Contrast audit of included transcriptomic accessions",
          (_p("literature/prisma/transcriptomics_contrast_audit.csv"),)),
    # New numbers, for files the text previously cited only by repository path.
    Table(11, "transcriptomic_screening", "Transcriptomic screening record",
          (_p("literature/prisma/transcriptomics_candidates.csv"),)),
    Table(12, "superseded_studies", "Study units superseded by a finer-grained re-analysis",
          (_p("conf/analysis/superseded_studies.csv"),)),
    Table(13, "phenotype_scope", "Phenotype scope decisions",
          (_p("conf/analysis/phenotype_scope.csv"),)),
    Table(14, "metabolite_id_map", "Metabolite name to ChEBI assignments",
          (_p("conf/analysis/metabolite_id_map.csv"),)),
    Table(15, "peptide_gene_map", "Peptide to precursor protein assignments",
          (_p("conf/analysis/peptide_gene_map.csv"),)),
    Table(16, "sni_tissue_annotations", "SNI tissue annotations",
          (_p("data/interim/sni_tissue_annotations.csv"),)),
    Table(17, "depth_sensitivity", "Subsampling to equal depth, per seed",
          (_p("results/meta/transcriptomics/depth_sensitivity.csv"),
           _p("results/meta/transcriptomics/depth_sensitivity_shared.csv"),
           _p("results/meta/transcriptomics/depth_sensitivity_strata.csv")),
          ("all_features", "shared_orthologs", "by_stratum")),
    # Table 5 of the Journal of Pain version. Brain Communications allows eight
    # display items, and Fig. 5A already shows the same ranking, so the table
    # moved to the supplementary (2026-09-21): an excerpt is printed there, and
    # this workbook carries the whole ranking.
    Table(18, "gwas_replication_ranking", "Lead-variant genes by replication count",
          (_p("results/meta/genomics/gwas_gene_meta_replication.csv"),)),
)


def _style(ws, n_cols: int) -> None:
    """Arial throughout, a bold frozen header, and columns wide enough to read."""
    for row in ws.iter_rows():
        for cell in row:
            cell.font = Font(name=FONT, size=10)
    for cell in ws[1]:
        cell.font = Font(name=FONT, size=10, bold=True)
        cell.alignment = Alignment(vertical="top", wrap_text=True)
    ws.freeze_panes = "A2"
    for i in range(1, n_cols + 1):
        letter = get_column_letter(i)
        widest = max((len(str(c.value)) for c in ws[letter] if c.value is not None),
                     default=10)
        ws.column_dimensions[letter].width = min(max(widest + 2, 10), MAX_COL_WIDTH)


def build(table: Table) -> Path:
    missing = [s for s in table.sources if not s.exists()]
    if missing:
        raise SystemExit(
            f"Table S{table.number}: missing source(s): "
            + ", ".join(str(m.relative_to(REPO_ROOT)) for m in missing))

    out = OUT_DIR / f"Table_S{table.number}_{table.slug}.xlsx"
    names = table.sheet_names()
    if len(names) != len(table.sources):
        raise SystemExit(f"Table S{table.number}: {len(names)} sheet names for "
                         f"{len(table.sources)} sources")

    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for src, sheet in zip(table.sources, names):
            df = pd.read_csv(src, dtype=str, keep_default_na=False)
            df.to_excel(writer, sheet_name=sheet[:31], index=False)
            _style(writer.sheets[sheet[:31]], df.shape[1])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", type=int, action="append",
                    help="build just this table number; repeatable")
    args = ap.parse_args()

    wanted = [t for t in TABLES if not args.only or t.number in args.only]
    if not wanted:
        raise SystemExit(f"no table matches --only {args.only}")

    total = 0
    for table in wanted:
        out = build(table)
        size = out.stat().st_size
        total += size
        print(f"  S{table.number:<2} {out.name:<46} {size/1e6:6.2f} MB")
    print(f"{len(wanted)} workbook(s), {total/1e6:.1f} MB -> "
          f"{OUT_DIR.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
