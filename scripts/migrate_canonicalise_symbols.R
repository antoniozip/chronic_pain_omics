#!/usr/bin/env Rscript
#' One-off migration: canonicalise gene symbols in existing normalized tables.
#'
#' 05b now canonicalises symbols as it writes (R/symbol_utils.R), but it cannot
#' be re-run on a workstation missing the Bioconductor annotation packages: the
#' lookup failure is silent and would replace gene symbols with probe ids for
#' every array study (see CLAUDE.md, Environment Caveats). This script applies
#' the identical repair to the tables already on disk, using the same shared
#' function, so the correction can land without re-deriving the annotation.
#'
#' It rewrites, in results/per_study/{modality}/:
#'   {accession}_effects_normalized.csv   — repaired symbols, merged duplicates
#'   effects_normalized_all.csv           — rebuilt from the repaired tables
#'
#' results/ is gitignored, so there is no version-control safety net: originals
#' are copied to a timestamped backup directory before anything is written, and
#' every substitution is recorded in a manifest.
#'
#' The repair is idempotent — a second run reports no changes.
#'
#' Usage:
#'   Rscript scripts/migrate_canonicalise_symbols.R --modality transcriptomics
#'   Rscript scripts/migrate_canonicalise_symbols.R --modality transcriptomics --apply
#'
#' Without --apply nothing is written; the script reports what it would do.

suppressPackageStartupMessages({
  library(dplyr)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))

source(file.path(REPO_ROOT, "R", "symbol_utils.R"))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
MODALITY <- {
  i <- match("--modality", args)
  if (!is.na(i) && length(args) > i) args[i + 1] else "transcriptomics"
}
APPLY <- "--apply" %in% args

IN_DIR <- file.path(REPO_ROOT, "results", "per_study", MODALITY)
if (!dir.exists(IN_DIR)) stop("No such directory: ", IN_DIR)

COMBINED <- file.path(IN_DIR, "effects_normalized_all.csv")
STAMP    <- format(Sys.time(), "%Y%m%d_%H%M%S")
BACKUP   <- file.path(IN_DIR, paste0(".backup_pre_canonicalise_", STAMP))
MANIFEST <- file.path(REPO_ROOT, "results", "per_study",
                      sprintf("canonicalisation_%s_%s.csv", MODALITY, STAMP))

message(sprintf("== Canonicalise symbols: %s ==", MODALITY))
message(if (APPLY) "Mode: APPLY (files will be rewritten)" else
        "Mode: dry run (nothing written; pass --apply to write)")

files <- list.files(IN_DIR, pattern = "_effects_normalized\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No normalized effect tables in ", IN_DIR)
message(sprintf("Found %d normalized tables", length(files)))

if (APPLY) dir.create(BACKUP, recursive = TRUE, showWarnings = FALSE)

all_changes <- list()
n_touched <- 0L
n_merged_total <- 0L

for (f in files) {
  accession <- sub("_effects_normalized\\.csv$", "", basename(f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$feature_id <- as.character(df$feature_id)
  if ("original_feature_id" %in% colnames(df)) {
    df$original_feature_id <- as.character(df$original_feature_id)
  }
  before <- nrow(df)

  res <- canonicalize_study_symbols(df, accession, verbose = TRUE)
  if (nrow(res$changes) == 0) next

  n_touched <- n_touched + 1L
  merged <- before - nrow(res$df)
  n_merged_total <- n_merged_total + merged
  all_changes[[accession]] <- res$changes

  if (APPLY) {
    file.copy(f, file.path(BACKUP, basename(f)), overwrite = FALSE)
    write.csv(res$df, f, row.names = FALSE)
  }
}

message(sprintf("\n%d of %d studies affected, %d row(s) merged",
                n_touched, length(files), n_merged_total))

if (length(all_changes) == 0) {
  message("Nothing to do — symbols are already canonical.")
  quit(save = "no", status = 0)
}

changes <- bind_rows(all_changes)
if (APPLY) {
  write.csv(changes, MANIFEST, row.names = FALSE)
  message(sprintf("Manifest → %s", MANIFEST))
  message(sprintf("Backup   → %s", BACKUP))
} else {
  message("\nSubstitutions that would be applied:")
  print(as.data.frame(changes %>% count(from, to, name = "studies") %>% arrange(from)))
}

# ---------------------------------------------------------------------------
# Rebuild the combined table from the repaired per-study tables
# ---------------------------------------------------------------------------
# 05b builds this by binding the studies it processed, which excludes the
# superseded ones; rebuild it the same way rather than patching the old file,
# so the combined table cannot disagree with its parts.

superseded <- superseded_studies(
  file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv"))
accs  <- sub("_effects_normalized\\.csv$", "", basename(files))
keep  <- files[!accs %in% superseded]
for (acc in accs[accs %in% superseded]) {
  message(sprintf("Combined table excludes %s: superseded", acc))
}

if (APPLY) {
  combined <- bind_rows(lapply(keep, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$feature_id <- as.character(d$feature_id)
    if ("original_feature_id" %in% colnames(d)) {
      d$original_feature_id <- as.character(d$original_feature_id)
    }
    d
  }))
  if (file.exists(COMBINED)) {
    file.copy(COMBINED, file.path(BACKUP, basename(COMBINED)), overwrite = FALSE)
  }
  write.csv(combined, COMBINED, row.names = FALSE)
  message(sprintf("Combined table rebuilt: %d rows from %d studies → %s",
                  nrow(combined), length(keep), basename(COMBINED)))
} else {
  message(sprintf("Would rebuild %s from %d studies", basename(COMBINED), length(keep)))
}

message("\nDownstream tables are now stale: re-run 06, 06b, 06g and 07.")
