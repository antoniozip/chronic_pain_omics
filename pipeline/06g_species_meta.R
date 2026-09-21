#!/usr/bin/env Rscript
#' Species-stratified meta-analysis.
#'
#' Partitions the per-study effect tables by species and runs a separate
#' random-effects meta-analysis within each. This is what makes a genuine
#' cross-species comparison possible: the whole-cohort pooled table
#' (06_meta_analysis.R) carries one estimate per feature over every study
#' regardless of species, so mouse and rat studies of the same gene symbol
#' collapse into a single number and cannot be told apart downstream.
#'
#' BH correction is applied within each species, over the features seen in at
#' least meta.min_studies studies of that species (conf/analysis/default.yaml);
#' features below that stay in the output with padj = NA. See adjust_pooled_fdr
#' in R/meta_utils.R.
#'
#' Studies contributing more than one species are excluded: the per-study
#' tables carry no sample-level species labels, so their rows cannot be
#' attributed. They are named in the log rather than dropped silently.
#'
#' Outputs:
#'   results/meta/transcriptomics/by_species/
#'     {species_slug}_pooled.csv   — pooled estimates within one species
#'     species_assignment.csv      — study → species, including exclusions
#'
#' Usage:
#'   Rscript pipeline/06g_species_meta.R --modality transcriptomics
#'   Rscript pipeline/06g_species_meta.R --modality transcriptomics --dry-run

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(jsonlite)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))

source(file.path(REPO_ROOT, "R", "meta_utils.R"))

args     <- commandArgs(trailingOnly = TRUE)
MODALITY <- {
  i <- match("--modality", args)
  if (!is.na(i) && length(args) > i) args[i + 1] else "transcriptomics"
}
DRY_RUN <- "--dry-run" %in% args

#' A sensitivity run needs to drop studies without pretending they are
#' superseded: the registry records a claim about sample independence, and
#' borrowing it to answer "what if this compartment were absent" would leave a
#' false statement in config. --exclude-studies takes its own list and
#' --out-dir keeps the result away from the published pools.
EXCLUDE_FILE <- {
  i <- match("--exclude-studies", args)
  if (!is.na(i) && length(args) > i) args[i + 1] else NA_character_
}
OUT_DIR_ARG <- {
  i <- match("--out-dir", args)
  if (!is.na(i) && length(args) > i) args[i + 1] else NA_character_
}

#' Study ids named in a one-column exclusion file, or character(0).
read_excluded_studies <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(character(0))
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (!"study_id" %in% colnames(d)) {
    stop("exclusion file needs a study_id column: ", path)
  }
  unique(trimws(as.character(d$study_id)))
}
EXCLUDED_STUDIES <- read_excluded_studies(EXCLUDE_FILE)

PER_STUDY_DIR <- file.path(REPO_ROOT, "results", "per_study", MODALITY)
MANIFEST      <- file.path(REPO_ROOT, "data", "interim", MODALITY, "harmonized.jsonl")
OUT_DIR       <- if (!is.na(OUT_DIR_ARG)) OUT_DIR_ARG else
  file.path(REPO_ROOT, "results", "meta", MODALITY, "by_species")

METHOD <- "REML"
MIN_STUDIES <- tryCatch(
  yaml::read_yaml(file.path(REPO_ROOT, "conf", "analysis", "default.yaml"))$meta$min_studies,
  error = function(e) 2L)
if (is.null(MIN_STUDIES)) MIN_STUDIES <- 2L

SPECIES <- c("Mus musculus", "Rattus norvegicus", "Homo sapiens")

species_slug <- function(x) gsub(" ", "_", x)


# ---------------------------------------------------------------------------
# Study → species assignment
# ---------------------------------------------------------------------------

#' Map each study id to its single species, or NA when it has none or several.
#'
#' Derived study ids (GSE241361_DRG, GSE241361_Spinal_cord) do not appear in
#' the manifest, which records only the parent accession. Their species is
#' inherited by stripping trailing `_<suffix>` components until an accession
#' matches, so a per-tissue split is not silently dropped from every pool.
build_species_map <- function(manifest_path, study_ids) {
  records <- lapply(readLines(manifest_path, warn = FALSE), function(l) {
    if (nzchar(trimws(l))) jsonlite::fromJSON(l) else NULL
  })
  records <- Filter(Negate(is.null), records)

  by_acc <- list()
  for (r in records) {
    sp <- r$species_canonical
    if (is.null(sp)) sp <- character(0)
    if (isTRUE(r$has_human) && !("Homo sapiens" %in% sp)) sp <- c(sp, "Homo sapiens")
    by_acc[[r$accession]] <- unique(as.character(sp))
  }

  resolve <- function(sid) {
    candidate <- sid
    while (nzchar(candidate)) {
      if (!is.null(by_acc[[candidate]])) return(by_acc[[candidate]])
      parent <- sub("_[^_]+$", "", candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
    character(0)
  }

  out <- lapply(study_ids, resolve)
  names(out) <- study_ids
  out
}


load_per_study_effects <- function(per_study_dir) {
  # Tiering must match 06's loader exactly: gene-keyed, then normalized, then
  # raw, one file per accession. If the two disagree the species pools are
  # computed from different tables than the whole-cohort pool.
  # `recursive` on every glob, because proteomics per-study tables live one
  # directory deeper (results/per_study/proteomics/<ACC>/<ACC>_effects.csv).
  gene_files <- list.files(per_study_dir, pattern = "_effects_gene\\.csv$",
                           full.names = TRUE, recursive = TRUE)
  norm_files <- list.files(per_study_dir, pattern = "_effects_normalized\\.csv$",
                           full.names = TRUE, recursive = TRUE)
  raw_files <- list.files(per_study_dir, pattern = "_effects\\.csv$", full.names = TRUE, recursive = TRUE)
  raw_files <- raw_files[!grepl("_normalized", raw_files)]

  gene_acc <- sub("_effects_gene\\.csv$", "", basename(gene_files))
  norm_acc <- sub("_effects_normalized\\.csv$", "", basename(norm_files))
  raw_acc  <- sub("_effects\\.csv$", "", basename(raw_files))
  norm_files <- norm_files[!norm_acc %in% gene_acc]
  norm_acc   <- norm_acc[!norm_acc %in% gene_acc]
  only_raw <- raw_files[!raw_acc %in% c(gene_acc, norm_acc)]

  files_to_load <- c(gene_files, norm_files, only_raw)
  # Same exclusion the whole-cohort meta applies: a finer-grained re-analysis
  # of the same samples supersedes the combined study.
  superseded <- superseded_studies(
    file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv"))
  if (length(superseded) > 0) {
    accs <- sub("_effects(_normalized|_gene)?\\.csv$", "", basename(files_to_load))
    for (acc in accs[accs %in% superseded]) {
      message(sprintf("  excluding %s: superseded", acc))
    }
    files_to_load <- files_to_load[!accs %in% superseded]
  }
  if (length(EXCLUDED_STUDIES) > 0) {
    accs <- sub("_effects(_normalized|_gene)?\\.csv$", "", basename(files_to_load))
    for (acc in accs[accs %in% EXCLUDED_STUDIES]) {
      message(sprintf("  excluding %s: named in --exclude-studies", acc))
    }
    files_to_load <- files_to_load[!accs %in% EXCLUDED_STUDIES]
  }
  if (length(files_to_load) == 0) stop("No per-study effect files in ", per_study_dir)

  bind_rows(lapply(files_to_load, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$feature_id <- as.character(d$feature_id)
    if ("original_feature_id" %in% colnames(d))
      d$original_feature_id <- as.character(d$original_feature_id)
    d
  }))
}


meta_one_feature <- function(sub, method = METHOD) {
  sub <- usable_se(sub)
  if (nrow(sub) == 0) return(NULL)
  study_ids <- paste(unique(sub$study_id), collapse = ";")
  if (nrow(sub) == 1) {
    r <- sub[1, ]
    return(data.frame(
      feature_id = r$feature_id, k = 1L, study_ids = study_ids,
      yi = r$effect_size, se = r$se,
      ci_lb = r$effect_size - Z_975 * r$se,
      ci_ub = r$effect_size + Z_975 * r$se,
      pval = r$pval, I2 = NA_real_, tau2 = NA_real_,
      stringsAsFactors = FALSE))
  }
  res <- tryCatch(
    rma(yi = effect_size, sei = se, data = sub, method = method, slab = sub$study_id),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  data.frame(
    feature_id = sub$feature_id[1], k = res$k, study_ids = study_ids,
    yi = res$b[1], se = res$se,
    ci_lb = res$ci.lb, ci_ub = res$ci.ub,
    pval = res$pval, I2 = res$I2, tau2 = res$tau2,
    stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

message(sprintf("[%s] BH correction family: features in >= %d studies of one species",
                MODALITY, MIN_STUDIES))

df <- load_per_study_effects(PER_STUDY_DIR)
df <- usable_se(df)
study_ids <- sort(unique(df$study_id))
message(sprintf("[%s] %d rows, %d studies, %d features",
                MODALITY, nrow(df), length(study_ids), length(unique(df$feature_id))))

species_map <- build_species_map(MANIFEST, study_ids)

# Curated per-study overrides. The manifest records what GEO says about a
# series; this records what was actually analysed, where the two differ. Only
# needed when a study's analysed subset is single-species while its GEO record
# is not -- otherwise the manifest governs. Mirrors superseded_studies.csv:
# an explicit, reasoned exception rather than editing the manifest, so the
# manifest keeps describing the source truthfully.
ovr_path <- file.path(REPO_ROOT, "conf", "analysis", "species_overrides.csv")
if (file.exists(ovr_path)) {
  ovr <- utils::read.csv(ovr_path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ovr))) {
    sid <- ovr$study_id[i]
    if (sid %in% names(species_map)) {
      message(sprintf("  species override: %s -> %s (conf/analysis/species_overrides.csv)",
                      sid, ovr$species[i]))
      species_map[[sid]] <- ovr$species[i]
    }
  }
}

n_species   <- vapply(species_map, length, integer(1))

assignment <- data.frame(
  study_id = study_ids,
  species  = vapply(species_map, function(s) paste(s, collapse = ";"), character(1)),
  status   = ifelse(n_species == 1, "assigned",
                    ifelse(n_species == 0, "excluded_no_species", "excluded_multi_species")),
  stringsAsFactors = FALSE
)
for (st in c("excluded_no_species", "excluded_multi_species")) {
  bad <- assignment[assignment$status == st, ]
  if (nrow(bad) > 0) {
    message(sprintf("  %s: %s", st, paste(bad$study_id, collapse = ", ")))
  }
}

df$species <- vapply(species_map[df$study_id],
                     function(s) if (length(s) == 1) s else NA_character_,
                     character(1))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(assignment, file.path(OUT_DIR, "species_assignment.csv"), row.names = FALSE)

for (sp in SPECIES) {
  sub_df <- df[!is.na(df$species) & df$species == sp, ]
  n_studies <- length(unique(sub_df$study_id))
  n_feats   <- length(unique(sub_df$feature_id))
  message(sprintf("\n=== %s: %d studies, %d features ===", sp, n_studies, n_feats))
  if (nrow(sub_df) == 0) next

  if (DRY_RUN) next

  feats    <- split(sub_df, sub_df$feature_id)
  res_list <- Filter(Negate(is.null), lapply(feats, meta_one_feature))
  pooled   <- bind_rows(res_list)
  if (nrow(pooled) == 0) next

  pooled <- pooled %>% arrange(pval) %>% mutate(species = sp)
  pooled <- adjust_pooled_fdr(pooled, min_k = MIN_STUDIES)

  out_csv <- file.path(OUT_DIR, sprintf("%s_pooled.csv", species_slug(sp)))
  write.csv(pooled, out_csv, row.names = FALSE)

  fam <- pooled[!is.na(pooled$padj), ]
  message(sprintf("  → %d features (%d in BH family, k>=%d)%s",
                  nrow(pooled), nrow(fam), MIN_STUDIES,
                  if (nrow(fam) > 0)
                    sprintf(", top hit: %s (padj=%.2e)", fam$feature_id[1], fam$padj[1])
                  else ""))
}

message("\nStep 06g complete.")
