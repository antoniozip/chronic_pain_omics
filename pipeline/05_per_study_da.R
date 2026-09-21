#!/usr/bin/env Rscript
#' Per-study differential analysis pipeline step.
#'
#' Reads harmonized manifests (data/interim/{modality}/harmonized.jsonl) and
#' the modality's own screening record (literature/prisma/{modality}_candidates.csv
#' or {modality}_manual_review.csv), then runs DA for each included study using
#' the appropriate modality-specific analysis module.
#'
#' The screening record must be keyed by dataset accession. It is deliberately
#' not literature/prisma/screening.csv, which is keyed by PubMed id — see
#' plan/2026-08-28-step05-screening-mismatch.md.
#'
#' Output: results/per_study/{modality}/{accession}_effects.csv
#'
#' Usage:
#'   Rscript pipeline/05_per_study_da.R --modality transcriptomics
#'   Rscript pipeline/05_per_study_da.R --modality all
#'   Rscript pipeline/05_per_study_da.R --modality transcriptomics --dry-run

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(jsonlite)
  library(dplyr)
})

# Source R modules relative to repo root
script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
REPO_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
R_DIR     <- file.path(REPO_ROOT, "R")

source(file.path(R_DIR, "da_utils.R"))
source(file.path(R_DIR, "transcriptomics_da.R"))
source(file.path(R_DIR, "genomics_da.R"))
source(file.path(R_DIR, "proteomics_da.R"))
source(file.path(R_DIR, "metabolomics_da.R"))
source(file.path(R_DIR, "single_cell_da.R"))

ALL_MODALITIES <- c("genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics")

# `single_cell` is analysable but is NOT in ALL_MODALITIES, so `--modality all`
# cannot reach it. That is the mechanism keeping the pseudobulked single-cell
# arm out of the bulk transcriptomic pool: merging the two requires someone to
# type the name, and there is no path by which it happens as a side effect.
# See plan/2026-09-08-single-cell-arm.md for why the arm is a stratum of its own.
EXTRA_MODALITIES <- c("single_cell")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

option_list <- list(
  make_option("--modality",       type = "character", default = "transcriptomics",
              help = "Modality or 'all' [default: transcriptomics]"),
  make_option("--analysis-config", type = "character",
              default = file.path(REPO_ROOT, "conf", "analysis", "default.yaml")),
  make_option("--ingest-config",   type = "character",
              default = file.path(REPO_ROOT, "conf", "ingest", "default.yaml")),
  make_option("--interim-dir",     type = "character",
              default = file.path(REPO_ROOT, "data", "interim")),
  make_option("--prisma-dir",      type = "character",
              default = file.path(REPO_ROOT, "literature", "prisma")),
  make_option("--out-dir",         type = "character",
              default = file.path(REPO_ROOT, "results", "per_study")),
  make_option("--cache-dir",       type = "character",
              default = file.path(REPO_ROOT, "data", "raw", "geo_cache")),
  make_option("--dry-run",         action = "store_true", default = FALSE),
  make_option("--skip-screening",  action = "store_true", default = FALSE,
              help = "Ignore screening decisions and process all harmonized records"),
  make_option("--accessions",      type = "character", default = "",
              help = paste("Comma-separated accessions to process, intersected with",
                           "the screening selection. For re-running specific studies",
                           "without recomputing the whole arm."))
)

opt <- parse_args(OptionParser(option_list = option_list))
cfg <- read_yaml(opt[["analysis-config"]])
modalities <- if (opt$modality == "all") ALL_MODALITIES else opt$modality
unknown <- setdiff(modalities, c(ALL_MODALITIES, EXTRA_MODALITIES))
if (length(unknown) > 0) {
  stop(sprintf("unknown modality: %s", paste(unknown, collapse = ", ")))
}


# ---------------------------------------------------------------------------
# Dispatch per modality
# ---------------------------------------------------------------------------

run_modality <- function(modality, cfg, opt) {
  message("\n========== ", toupper(modality), " ==========")

  harmonized   <- load_harmonized(opt[["interim-dir"]], modality)
  if (nrow(harmonized) == 0) return(invisible(NULL))

  if (!isTRUE(opt[["skip-screening"]])) {
    included_ids <- load_included_accessions(opt[["prisma-dir"]], modality)
    # Refuse to narrow the corpus with a list that cannot select from it.
    assert_screening_identifier_space(included_ids, harmonized$accession, modality)
    if (!is.null(included_ids) && length(included_ids) > 0) {
      n_before   <- nrow(harmonized)
      keep       <- accession_screened_in(harmonized$accession, included_ids)
      inherited  <- harmonized$accession[keep &
                                         !(harmonized$accession %in% included_ids)]
      harmonized <- harmonized[keep, ]
      message(sprintf("[%s] Screening filter: %d of %d harmonized records retained",
                      modality, nrow(harmonized), n_before))
      if (length(inherited) > 0) {
        message(sprintf("[%s] Verdict inherited from parent accession: %s",
                        modality, paste(sort(inherited), collapse = ", ")))
      }
    }
  }

  requested <- trimws(strsplit(opt[["accessions"]], ",")[[1]])
  requested <- requested[nzchar(requested)]
  if (length(requested) > 0) {
    absent <- setdiff(requested, harmonized$accession)
    if (length(absent) > 0) {
      # Silently processing fewer studies than named is the failure mode this
      # whole script has already been bitten by twice. Name them.
      message(sprintf("[%s] Not selected, so not processed: %s",
                      modality, paste(absent, collapse = ", ")))
    }
    harmonized <- harmonized[harmonized$accession %in% requested, ]
    if (nrow(harmonized) == 0) {
      stop(sprintf("[%s] none of the requested accessions are in the selection",
                   modality))
    }
  }

  message(sprintf("[%s] %d studies to process", modality, nrow(harmonized)))
  if (opt[["dry-run"]]) {
    message(sprintf("[dry-run][%s] Would process: %s",
                    modality, paste(harmonized$accession, collapse = ", ")))
    return(invisible(NULL))
  }

  results <- vector("list", nrow(harmonized))
  for (i in seq_len(nrow(harmonized))) {
    row <- harmonized[i, ]
    acc <- row$accession

    effect_df <- tryCatch(
      switch(modality,
        transcriptomics = run_transcriptomics_da(
          acc, row$platform_canonical, cfg, cache_dir = opt[["cache-dir"]]
        ),
        genomics = run_genomics_da(acc, row$effect_size_unit, cfg),
        proteomics = {
          message("  [proteomics] Raw data download for PRIDE requires manual setup — skipping ", acc)
          NULL
        },
        metabolomics = ,
        lipidomics   = {
          message("  [", modality, "] Raw data download for MetaboLights requires manual setup — skipping ", acc)
          NULL
        },
        single_cell = run_single_cell_da(acc, cfg, interim_dir = opt[["interim-dir"]])
      ),
      error = function(e) {
        message("  ERROR for ", acc, ": ", e$message)
        NULL
      }
    )

    if (!is.null(effect_df) && nrow(effect_df) > 0) {
      save_effects(effect_df, opt[["out-dir"]], modality, acc)
      results[[i]] <- effect_df
    }
  }

  completed <- Filter(Negate(is.null), results)
  message(sprintf("[%s] Completed DA for %d / %d studies",
                  modality, length(completed), nrow(harmonized)))
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

for (mod in modalities) {
  run_modality(mod, cfg, opt)
}

message("\nStep 05 complete. Effect files written to: ", opt[["out-dir"]])
