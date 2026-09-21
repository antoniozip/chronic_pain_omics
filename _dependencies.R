#' Dependencies renv cannot discover by static analysis.
#'
#' This file is never executed. It exists so `renv::snapshot()`, which scans
#' the project for `library()`/`require()` calls, records packages that the
#' pipeline loads dynamically and would otherwise omit from renv.lock.
#'
#' pipeline/05b_normalize_ids.R selects a Bioconductor annotation package at
#' run time from the GEO platform id:
#'
#'     id_map <- map_via_annotation_db(df$feature_id, GPL_TO_ANNO_PKG[[gpl]])
#'
#' The package name is a string looked up in a table, so nothing in the source
#' names it literally. rat2302.db and rgu34a.db were missing from renv.lock for
#' exactly this reason; without them 05b does not fail, it silently falls
#' through to "keeping native IDs" and nine studies degrade from gene symbols
#' to raw probe ids. Add an entry here whenever GPL_TO_ANNO_PKG gains one.

#' Every entry of GPL_TO_ANNO_PKG must appear below. Only 5 of the 11 did until
#' 2026-08-15, so `renv::restore()` could not install the other 6 packages and a
#' container built from the lockfile would have reproduced the silent
#' probe-id fallback this file exists to prevent. `containers/verify_environment.R`
#' now checks the two lists against each other so they cannot drift again.

# Affymetrix rat arrays
library(rgu34a.db)                   # GPL85    - Rat Genome U34A
library(rgu34b.db)                   # GPL86    - Rat Genome U34B
library(rgu34c.db)                   # GPL87    - Rat Genome U34C
library(rae230a.db)                  # GPL341   - Rat Expression 230A
library(rat2302.db)                  # GPL1355  - Rat Genome 230 2.0

# Affymetrix mouse arrays
library(mouse4302.db)                # GPL1261  - Mouse Genome 430 2.0
library(mogene10sttranscriptcluster.db)  # GPL6246/GPL6247 - Mouse Gene 1.0 ST
library(mogene11sttranscriptcluster.db)  # GPL16570        - Mouse Gene 1.1 ST
library(mta10transcriptcluster.db)   # mouse MTA-1.0 (also loaded literally)

# Illumina BeadArrays
library(illuminaMousev1.db)          # GPL6194/GPL6543 - MouseRef-8 v1.1
library(illuminaMousev2.db)          # GPL6885/GPL6887 - MouseRef-8 v2.0
library(illuminaRatv1.db)            # GPL2050  - RatRef-12 v1.0

# Organism-level symbol mapping. Reached through
# `getExportedValue(pkg, pkg)` with the package name looked up from a species
# tag, so nothing in 05b names them literally -- the same blind spot as
# GPL_TO_ANNO_PKG. map_entrez() and map_ensembl_via_orgdb() both depend on
# these; without them GSE117526 keeps bare Entrez ids and the two human
# studies keep bare ENSG ids, silently, since both routes just return NULL.
library(org.Mm.eg.db)                # mouse gene symbol mapping
library(org.Rn.eg.db)                # rat gene symbol mapping
library(org.Hs.eg.db)                # human gene symbol mapping

#' biomaRt is the fallback when GPL_TO_ANNO_PKG has no entry for a platform.
#' Like the annotation packages it is reached through a runtime branch rather
#' than a literal library() call in 05b, so renv's scan misses it; without it
#' the fallback silently keeps native probe ids, the same failure mode the
#' annotation packages have.
library(biomaRt)
