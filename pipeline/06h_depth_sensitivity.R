#!/usr/bin/env Rscript
#' Is the mouse/rat significance gap an artefact of sampling depth?
#'
#' Mouse and rat contribute 26 studies each, yet the mouse pool yields almost
#' no significant feature (8 of 35,459) where rat yields 876 of 47,929. The
#' obvious suspect is depth: features in the mouse correction family are seen
#' in a median of 13 mouse studies against 3 for rat, and a random-effects
#' pool over more studies exposes more between-study heterogeneity.
#'
#' This script tests that explanation two ways:
#'
#'   1. Matched depth. Compare the significant fraction within the k=3 stratum
#'      of each species, where no subsampling is needed at all.
#'   2. Subsampling, applied symmetrically. Draw 3 studies at random for every
#'      feature with k>=3, re-pool, re-apply BH over the resulting family, and
#'      count; repeat over several seeds. Both species are subsampled, not only
#'      the deeper one: comparing a subsampled mouse against an observed rat
#'      would leave the rat side untouched by the very procedure whose effect
#'      is being measured. Rat is barely affected in practice, since its median
#'      k is already 3, which is itself the point.
#'
#' Species assignment is read from 06g's species_assignment.csv rather than
#' re-derived, so the two cannot disagree.
#'
#' Output:
#'   results/meta/{modality}/depth_sensitivity.csv       — per-seed counts
#'   results/meta/{modality}/depth_sensitivity_strata.csv — observed by k bin
#'
#' Usage:
#'   Rscript pipeline/06h_depth_sensitivity.R --modality transcriptomics
#'   Rscript pipeline/06h_depth_sensitivity.R --modality transcriptomics --seeds 5

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) > i) args[i + 1] else default
}
MODALITY <- arg_val("--modality", "transcriptomics")
SPECIES  <- c("Mus musculus", "Rattus norvegicus")
N_SEEDS  <- as.integer(arg_val("--seeds", "5"))
TARGET_K <- as.integer(arg_val("--target-k", "3"))

PER_STUDY_DIR <- file.path(REPO_ROOT, "results", "per_study", MODALITY)
BY_SPECIES    <- file.path(REPO_ROOT, "results", "meta", MODALITY, "by_species")
OUT_DIR       <- file.path(REPO_ROOT, "results", "meta", MODALITY)

METHOD <- "REML"
MIN_STUDIES <- tryCatch(
  yaml::read_yaml(file.path(REPO_ROOT, "conf", "analysis", "default.yaml"))$meta$min_studies,
  error = function(e) 3L)

assign_path <- file.path(BY_SPECIES, "species_assignment.csv")
if (!file.exists(assign_path)) {
  stop("Run pipeline/06g_species_meta.R first: no ", assign_path)
}
assignment <- read.csv(assign_path, stringsAsFactors = FALSE)
assignment <- assignment[assignment$status == "assigned", ]



# ---------------------------------------------------------------------------
# 1. Observed significant fraction by depth, both species
# ---------------------------------------------------------------------------

strata <- list()
for (sp in c("Mus musculus", "Rattus norvegicus")) {
  f <- file.path(BY_SPECIES, sprintf("%s_pooled.csv", gsub(" ", "_", sp)))
  d <- read.csv(f, stringsAsFactors = FALSE)
  fam <- d[!is.na(d$padj), ]
  bins <- cut(fam$k, c(2, 3, 4, 6, 10, 20, Inf),
              labels = c("3", "4", "5-6", "7-10", "11-20", "21+"))
  s <- data.frame(species = sp, k_bin = bins, sig = fam$padj < 0.05,
                  stringsAsFactors = FALSE) %>%
    group_by(species, k_bin) %>%
    summarise(n = n(), n_sig = sum(sig), .groups = "drop") %>%
    mutate(pct_sig = 100 * n_sig / n)
  strata[[sp]] <- s
  message(sprintf("  %s: family %d, significant %d (%.2f%%), median k %.0f",
                  sp, nrow(fam), sum(fam$padj < 0.05),
                  100 * mean(fam$padj < 0.05), median(fam$k)))
}
strata_df <- bind_rows(strata)
write.csv(strata_df, file.path(OUT_DIR, "depth_sensitivity_strata.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. Subsample mouse to rat-like depth
# ---------------------------------------------------------------------------

load_studies <- function(dir, ids) {
  files <- c(
    file.path(dir, paste0(ids, "_effects_normalized.csv")),
    file.path(dir, paste0(ids, "_effects.csv"))
  )
  files <- files[file.exists(files)]
  seen <- character(0)
  out <- list()
  for (f in files) {
    acc <- sub("_effects(_normalized)?\\.csv$", "", basename(f))
    if (acc %in% seen) next          # prefer the normalized table
    seen <- c(seen, acc)
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$feature_id <- as.character(d$feature_id)
    out[[acc]] <- d[, intersect(colnames(d),
                                c("feature_id", "effect_size", "se", "pval", "study_id"))]
  }
  bind_rows(out)
}

pool_rows <- function(sub) {
  if (nrow(sub) == 1) {
    return(c(yi = sub$effect_size[1], pval = sub$pval[1]))
  }
  res <- tryCatch(rma(yi = effect_size, sei = se, data = sub, method = METHOD),
                  error = function(e) NULL)
  if (is.null(res)) return(c(yi = NA_real_, pval = NA_real_))
  c(yi = res$b[1], pval = res$pval)
}

#' Subsample one species to TARGET_K studies per feature and re-pool.
#'
#' Run for both species, not only the deeper one: comparing a subsampled mouse
#' against an observed rat would leave the rat side untouched by the very
#' procedure whose effect is being measured.
subsample_species <- function(sp, restrict_to = NULL, label = "all") {
  studies <- assignment$study_id[assignment$species == sp]
  message(sprintf("\n== %s [%s]: %d studies ==", sp, label, length(studies)))

  df <- usable_se(load_studies(PER_STUDY_DIR, studies))
  counts <- table(df$feature_id)
  eligible <- names(counts)[counts >= MIN_STUDIES]
  if (!is.null(restrict_to)) eligible <- intersect(eligible, restrict_to)
  message(sprintf("  rows with usable SE: %d | features with k>=%d: %d",
                  nrow(df), MIN_STUDIES, length(eligible)))

  df <- df[df$feature_id %in% eligible, ]
  split_feats <- split(df, df$feature_id)
  n_deeper <- sum(vapply(split_feats, function(x) nrow(x) > TARGET_K, logical(1)))
  message(sprintf("  features deeper than k=%d (actually subsampled): %d (%.1f%%)",
                  TARGET_K, n_deeper, 100 * n_deeper / length(split_feats)))

  out <- list()
  for (seed in seq_len(N_SEEDS)) {
    set.seed(1000 + seed)
    t0 <- Sys.time()
    pv <- vapply(split_feats, function(sub) {
      idx <- if (nrow(sub) > TARGET_K) sample.int(nrow(sub), TARGET_K) else seq_len(nrow(sub))
      pool_rows(sub[idx, , drop = FALSE])["pval"]
    }, numeric(1))

    padj <- p.adjust(pv, method = "BH")
    n_sig <- sum(padj < 0.05, na.rm = TRUE)
    n_tot <- sum(!is.na(padj))
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    message(sprintf("  seed %d: %d of %d significant (%.3f%%) [%s min]",
                    seed, n_sig, n_tot, 100 * n_sig / n_tot, el))
    out[[seed]] <- data.frame(
      species = sp, feature_set = label, seed = seed, target_k = TARGET_K,
      n_features = n_tot, n_subsampled = n_deeper,
      n_sig = n_sig, pct_sig = 100 * n_sig / n_tot
    )
  }
  bind_rows(out)
}

res_df <- bind_rows(lapply(SPECIES, subsample_species))
write.csv(res_df, file.path(OUT_DIR, "depth_sensitivity.csv"), row.names = FALSE)

message(sprintf("\n== Subsampled to k=%d, mean over %d seeds ==", TARGET_K, N_SEEDS))
for (sp in SPECIES) {
  r <- res_df[res_df$species == sp, ]
  obs <- strata_df[strata_df$species == sp, ]
  obs_all <- sum(obs$n_sig) / sum(obs$n) * 100
  message(sprintf("  %-20s observed %.3f%%  ->  subsampled %.3f%% (range %.3f-%.3f)",
                  sp, obs_all, mean(r$pct_sig), min(r$pct_sig), max(r$pct_sig)))
}


# ---------------------------------------------------------------------------
# 3. Matched composition: features whose human ortholog is seen in both species
# ---------------------------------------------------------------------------
# The two correction families cover different feature sets, so part of the
# residual gap at common depth could be composition rather than precision.
# Restricting both species to the animal features whose human ortholog appears
# in both concordance tables removes that difference. BH is re-applied within
# the restricted set, which is a self-contained family.

CS_DIR <- file.path(REPO_ROOT, "results", "cross_species", MODALITY)
conc_files <- c("Mus musculus" = "concordance_musculus.csv",
                "Rattus norvegicus" = "concordance_norvegicus.csv")

if (!all(file.exists(file.path(CS_DIR, conc_files)))) {
  message("\nNo concordance tables (run pipeline/07_cross_species.py) — ",
          "skipping the shared-ortholog analysis.")
} else {
  conc <- lapply(conc_files, function(f)
    read.csv(file.path(CS_DIR, f), stringsAsFactors = FALSE))
  names(conc) <- names(conc_files)
  shared_human <- Reduce(intersect, lapply(conc, function(d) unique(d$human_feature_id)))
  message(sprintf("\nHuman orthologs seen in both species: %d", length(shared_human)))

  restrict <- lapply(conc, function(d)
    unique(d$animal_feature_id[d$human_feature_id %in% shared_human]))

  # Observed rates on the shared set, BH re-applied within it
  obs_shared <- list()
  for (sp in SPECIES) {
    f <- file.path(BY_SPECIES, sprintf("%s_pooled.csv", gsub(" ", "_", sp)))
    d <- read.csv(f, stringsAsFactors = FALSE)
    d <- d[d$feature_id %in% restrict[[sp]] & !is.na(d$padj), ]
    padj <- p.adjust(d$pval, method = "BH")
    obs_shared[[sp]] <- data.frame(
      species = sp, feature_set = "shared_orthologs", seed = NA_integer_,
      target_k = NA_integer_, n_features = nrow(d), n_subsampled = NA_integer_,
      n_sig = sum(padj < 0.05, na.rm = TRUE),
      pct_sig = 100 * mean(padj < 0.05, na.rm = TRUE)
    )
    message(sprintf("  %-20s observed on shared set: %d of %d (%.3f%%), median k %.0f",
                    sp, obs_shared[[sp]]$n_sig, nrow(d), obs_shared[[sp]]$pct_sig,
                    median(d$k)))
  }

  sub_shared <- bind_rows(lapply(SPECIES, function(sp)
    subsample_species(sp, restrict_to = restrict[[sp]], label = "shared_orthologs")))

  shared_df <- bind_rows(bind_rows(obs_shared), sub_shared)
  write.csv(shared_df, file.path(OUT_DIR, "depth_sensitivity_shared.csv"),
            row.names = FALSE)

  message("\n== Shared ortholog set, subsampled to k=3 ==")
  for (sp in SPECIES) {
    r <- sub_shared[sub_shared$species == sp, ]
    o <- obs_shared[[sp]]
    message(sprintf("  %-20s observed %.3f%%  ->  subsampled %.3f%% (range %.3f-%.3f)",
                    sp, o$pct_sig, mean(r$pct_sig), min(r$pct_sig), max(r$pct_sig)))
  }
}

message("\nStep 06h complete.")
