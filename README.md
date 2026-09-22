# chronic_pain_omics

Analysis code for **"Meta-analysis of chronic pain multi-omics datasets reveals disjoint evidence"** (A. G. Zippo, 2026): a systematic, PRISMA-style multi-omics meta-analysis of
**chronic pain**, integrating genomics, transcriptomics, proteomics and
metabolomics across human cohorts and rodent models.

This is a code-only snapshot. The manuscript sources are not included, so the
manuscript claim tests (`tests/test_manuscript_claims.py`) skip here, and
paths under `manuscript/` mentioned below are written by the pipeline rather
than shipped with it.

This repository holds the analysis code, the configuration and the screening
records. It is not a model-training project: it is
literature screening, data harmonisation across heterogeneous omics
repositories, statistical pooling, cross-species integration, and the
bookkeeping that makes those reproducible.

**Paper:** *Meta-analysis of chronic pain multi-omics datasets reveals disjoint evidence* (`manuscript/manuscript.pdf`).

## What it found

| | |
| --- | --- |
| Transcriptomic study units pooled | 93 |
| Features with a pooled estimate | 161,284 (93,861 in the correction family) |
| GWAS with harmonised summary statistics | 57 |
| Proteomic datasets screened / poolable | 43 / 3 |
| Metabolomic study units | 12 |
| Cross-species ortholog pairs | 33,029 |

Three results, all negative, and that is the point:

- **Pain model adds modest structure, no more.** 11,043 of 72,236 testable
  features differ significantly between pain models, but shuffling the model
  label across studies still makes 7.8% of features look model-dependent. The
  excess over that null, roughly 5,600 features, falls short of conventional
  significance (p = 0.066, 1,000 permutations).
- **Rodent-to-human agreement is at or below chance.** Direction concordance
  is 51.0% for mouse and 52.2% for rat, falling to 49.0% and 50.7% when the
  human pool is restricted to tissue-comparable studies. Two independent
  halves of one human condition, through the same pipeline, agree 55.2%, so
  the cross-species null measures the reproducibility floor of pooled bulk
  transcriptomes rather than a fact about species.
- **Genomic and transcriptomic evidence implicate disjoint genes.** None of the
  25 genome-wide-significant genes is significant in the human transcriptomic
  pool or is a cross-species concordant transcript, and none of the 15
  most-replicated GWAS pain genes approaches significance there (smallest
  adjusted p 0.27).

## Layout

```
conf/          YAML configuration: search queries, inclusion criteria, analysis knobs
literature/    PRISMA records: every screening decision with a reason code
pipeline/      Numbered orchestration scripts, 01 through 15
src/           Python package (ingestion, harmonisation, plotting)
R/             R helpers (per-study DA, meta-analysis, shared statistics)
scripts/       Manuscript artefacts and maintenance tools
tests/         pytest and testthat suites
manuscript/    LaTeX sources, figures, supplementary
results/       Analysis outputs (gitignored: 12 GB, regenerable)
```

## Running it

```bash
uv venv && source .venv/bin/activate && uv sync   # Python
R -e 'renv::restore()'                            # R
```

The pipeline runs in order; each step reads the previous step's outputs from
`results/`.

```bash
python pipeline/01_search_literature.py --query-config conf/search/chronic_pain.yaml
python pipeline/02_screen.py
python pipeline/03_ingest_omics.py --modality transcriptomics
python pipeline/04_harmonize.py --modality transcriptomics
Rscript pipeline/05_per_study_da.R --modality transcriptomics
Rscript pipeline/06_meta_analysis.R --modality transcriptomics
Rscript pipeline/06g_species_meta.R --modality transcriptomics
python pipeline/07_cross_species.py
python pipeline/08_figures.py
```

Transcriptomics is the slow arm: steps 05 and 06 take roughly two hours each.
`CLAUDE.md` documents the caveats that are easy to get wrong, and is worth
reading before changing any analysis step.

## The manuscript's numbers are tested

Every numeric claim in the paper is bound to the results file that produces it.

```bash
pytest tests/test_manuscript_claims.py          # 326 claims against results/
python scripts/mutation_test_claims.py          # prove each claim can fail
```

The registry pairs a regex locating the sentence with a function computing the
value from the tracked outputs. A rewritten sentence whose pattern no longer
matches **fails** rather than skipping, and a figure stated twice fails too,
because the same number written in two places is how a Discussion drifts away
from its Results. Three rounds of manual review had still left four stale
numbers behind; the registry caught all four in one run.

Claims are themselves mutation-tested: `mutation_test_claims.py` perturbs the
number each one reads and checks that the claim fails. A claim that survives
its own mutation reports green while checking nothing. **326 of 326 are
killed.**

## Building the paper

```bash
cd manuscript
pdflatex manuscript && bibtex manuscript && pdflatex manuscript && pdflatex manuscript
pdflatex supplementary && pdflatex supplementary   # after manuscript.aux exists
```

References are BibTeX (`manuscript/references.bib`) in Elsevier's
`elsarticle-num` style, which ships in TeX Live's `texlive-publishers`. Every
entry carries the DOI that was resolved against CrossRef: when the hand-written
bibliography was converted, ten of its 46 entries did not survive that check,
four of them naming papers that do not exist under the title given.

`python scripts/build_manuscript_docx.py` renders Word versions.

## Tests

```bash
pytest tests/                                         # 782
Rscript -e 'testthat::test_dir("tests/testthat")'     # 306
ruff check . && mypy src/
```

CI runs all three on every push, plus the claim registry against a cached copy
of `results/`. A cache miss fails the job by design: `results/` cannot be
rebuilt inside a job, and an unverified manuscript should not be able to report
green.

## Data availability

No primary data is generated here. Every dataset is retrieved by the pipeline
from public repositories: GEO, the GWAS Catalog, PRIDE, MetaboLights and
Metabolomics Workbench. Accessions for every included dataset are in
`manuscript/supplementary/`, and every inclusion and exclusion decision is
recorded with a reason code in `literature/prisma/`.

`results/` is gitignored as regenerable, with one deliberate exception:
`results/meta/transcriptomics/egger_test.csv`, which the manuscript's
publication-bias claims read.

## Citing and licence

If you use this code or the harmonised effect-size tables, please cite the
paper. Code is MIT; text and data are CC-BY 4.0.

Maintainer: Antonio G. Zippo, Institute of Neuroscience, Consiglio Nazionale
delle Ricerche, and Milan Center for Neuroscience (NeuroMI).
