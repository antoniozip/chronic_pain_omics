# Analysis container

Closes the two reproducibility gaps recorded in
`plan/2026-08-15-reproducibility-check.md`.

| gap | before | after |
|---|---|---|
| R environment **described but absent** | `renv.lock` pins R 4.5.3 / Bioconductor 3.22; the workstation runs R 4.3.3 with 56 of 145 locked packages missing, so `05b` and `05` cannot run | the base image *is* R 4.5 / Bioc 3.22, and every package `GPL_TO_ANNO_PKG` can select is installed |
| MAGMA **present but undescribed** | binary and 1000G EUR panel under a gitignored `tools/`, captured by no lockfile; a fresh clone had no reproducible way to obtain them | binary baked in against a recorded checksum; panel fetched by a script that verifies its own download |

## Build

```bash
docker build -f containers/Dockerfile -t cp-multiomics:latest .
```

The build context is ~400 MB: `.dockerignore` excludes `results/` (12.6 GB) and
`tools/` (4 GB), both of which are mounted at run time instead.

## Verify

This is the step that matters, and it is written to fail rather than reassure:

```bash
docker run --rm cp-multiomics:latest Rscript containers/verify_environment.R
```

It checks the running R and Bioconductor versions against `renv.lock`, parses
`GPL_TO_ANNO_PKG` out of `pipeline/05b_normalize_ids.R` and requires every
package it names to be loadable, requires `_dependencies.R` to cover the same
set, and requires `magma` on `PATH`. Run on the host that motivated all this,
it exits 1 with nine absent annotation packages — which is the point.

**Why annotation packages get their own check.** A missing one does not make
`05b` fail. It falls through to keeping the array's native probe identifiers,
nine studies silently degrade from gene symbols to probe ids, and every
downstream number moves for a reason nobody can see. The same silence is why
`_dependencies.R` is checked against `GPL_TO_ANNO_PKG`: renv discovers
dependencies by static analysis and `05b` selects its package from a lookup
table, so the shim is the only thing that puts them in the lockfile — and it
listed 5 of 11 until 2026-08-15.

## Reference data

The 1000 Genomes EUR panel is ~4 GB and is deliberately not in the image:

```bash
containers/fetch_reference_data.sh              # download and unpack
containers/fetch_reference_data.sh --verify-only
```

Checksums are measured from the artefacts that produced the published results,
so a changed upstream archive fails loudly instead of quietly altering the
gene-based analysis.

## Run the pipeline

Mount the large directories rather than copying them in:

```bash
docker run --rm -it \
  -v "$PWD/data:/work/data" \
  -v "$PWD/results:/work/results" \
  -v "$PWD/tools:/work/tools" \
  cp-multiomics:latest bash
```

Then, inside:

```bash
Rscript pipeline/05b_normalize_ids.R          # the step the host cannot run
Rscript pipeline/05_per_study_da.R --modality transcriptomics
Rscript pipeline/06_meta_analysis.R --modality transcriptomics
uv run python pipeline/07_cross_species.py --modality transcriptomics
```

**Spot-check `05b` before trusting a full re-run.** Confirm a repaired study's
`feature_id` column holds gene symbols rather than probe ids; the failure this
container exists to prevent is invisible in the exit code.

## What the container does not fix

It pins the software, not the inputs. GEO, PRIDE and MetaboLights are queried
live by steps 01–05, and repositories add and revise datasets, so a rerun a
year from now may retrieve a different corpus. The per-study effect-size tables
are the stable entry point; the manifests under `literature/prisma/` record
which accessions the published analysis used.
