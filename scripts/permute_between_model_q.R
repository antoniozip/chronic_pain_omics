#!/usr/bin/env Rscript
#' Permutation null for the between-model Q test.
#'
#' 06b reports that 14,267 of 67,352 features (21.2%) show significant
#' between-model differences, and the paper reads that as pain-model
#' specificity. The reviewer objection -- and, after the 2026-09-03 positive
#' control showed effect estimates do not reproduce between two halves of a
#' single condition, the live concern -- is that Q_between is significant
#' because *studies* differ, not because *models* do. With non-reproducing
#' estimates, any partition of the studies would separate them.
#'
#' This tests it directly. The model label is permuted across study units,
#' keeping each study's effect vector intact: platform, tissue, sample size,
#' depth and small-study character are preserved exactly, and only the
#' assignment of studies to models is destroyed.
#'
#' Permutation is **within species**. Models nest cleanly inside species here
#' (all five rodent models are rodent, every human stratum is human), so a free
#' shuffle would put human studies into "CCI" and rodent studies into
#' "endometriosis_human". That destroys a constraint that is not under test and
#' would give an artificially easy null. Shuffling within species preserves the
#' species/model alignment and the size of every stratum, and asks only whether
#' *which* model a study belongs to carries information.
#'
#' Features are subsampled for cost: the full test is ~90 min over 67,352
#' features. The same subsample is used for the observed and every permuted
#' statistic, and BH correction is applied within the subsample each time, so
#' the counts are directly comparable.
#'
#' Usage:
#'   Rscript scripts/permute_between_model_q.R --features 3000 --perms 20

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args); if (!is.na(i) && length(args) > i) args[i + 1] else default
}
N_FEATURES <- as.integer(argval("--features", "3000"))
N_PERMS    <- as.integer(argval("--perms", "20"))
SEED       <- as.integer(argval("--seed", "1"))
#' The feature subsample is seeded separately from the permutations so a long
#' run can be split across processes: same features, disjoint label shuffles,
#' results poolable. Permutation seeds are SEED*1000 + p, so distinct --seed
#' values cannot collide for runs of fewer than 1000 permutations.
FEATURE_SEED <- as.integer(argval("--feature-seed", argval("--seed", "1")))
OUT        <- argval("--out", file.path(REPO_ROOT, "results", "meta",
                                        "transcriptomics", "stratified",
                                        "between_model_Q_permutation.csv"))
METHOD <- "REML"

# --- the stratum map, exactly as 06b builds it -----------------------------
src <- readLines(file.path(REPO_ROOT, "pipeline", "06b_stratified_meta.R"))
txt <- paste(src, collapse = "\n")
blk <- sub("^.*?STUDY_MODEL <- c\\(", "", txt)
blk <- sub("\\n\\).*$", "", blk)
pairs <- regmatches(blk, gregexpr('(GSE[A-Za-z0-9_]+)\\s*=\\s*"[^"]+"', blk))[[1]]
STUDY_MODEL <- setNames(sub('.*"([^"]+)".*', "\\1", pairs),
                        sub('^\\s*(GSE[A-Za-z0-9_]+)\\s*=.*', "\\1", pairs))
sf <- read.csv(file.path(REPO_ROOT, "conf", "analysis", "study_strata.csv"),
               stringsAsFactors = FALSE)
extra <- setNames(as.character(sf$stratum), sf$study_id)
STRATUM_MAP <- STUDY_MODEL
STRATUM_MAP[names(extra)] <- extra
SUBGROUPS <- unique(c("CCI", "SNI", "SNL", "CFA", "CIPN", "neuropathic_human",
                      "lbp_human", "nociplastic_human", "invitro_human",
                      unname(extra)))

# --- data ------------------------------------------------------------------
EFFECTS <- file.path(REPO_ROOT, "results", "per_study", "transcriptomics",
                     "effects_normalized_all.csv")
df <- read.csv(EFFECTS, stringsAsFactors = FALSE)
df <- usable_se(df)
sup <- superseded_studies(file.path(REPO_ROOT, "conf", "analysis",
                                    "superseded_studies.csv"))
df <- df[!(df$study_id %in% sup), ]
df$pain_model <- resolve_by_parent(df$study_id, STRATUM_MAP)
df <- df[df$pain_model %in% SUBGROUPS, ]

spa <- read.csv(file.path(REPO_ROOT, "results", "meta", "transcriptomics",
                          "by_species", "species_assignment.csv"),
                stringsAsFactors = FALSE)
SPECIES_OF <- setNames(spa$species, spa$study_id)

units <- unique(df$study_id)
unit_model   <- setNames(df$pain_model[match(units, df$study_id)], units)
unit_species <- SPECIES_OF[units]
keep <- !is.na(unit_species)
units <- units[keep]; unit_model <- unit_model[keep]; unit_species <- unit_species[keep]
df <- df[df$study_id %in% units, ]
message(sprintf("%d rows, %d study units, %d models, %d species",
                nrow(df), length(units), length(unique(unit_model)),
                length(unique(unit_species))))

# --- eligible features, subsampled ----------------------------------------
elig <- df %>%
  group_by(feature_id) %>%
  summarise(n_models = n_distinct(pain_model), k_total = n(), .groups = "drop") %>%
  filter(n_models >= 2, k_total >= 4)
message(sprintf("eligible features (>=2 models, k>=4): %d", nrow(elig)))
set.seed(FEATURE_SEED)
feats <- sample(elig$feature_id, min(N_FEATURES, nrow(elig)))
sub_df <- df[df$feature_id %in% feats, c("feature_id", "study_id",
                                         "effect_size", "se")]
by_feat <- split(sub_df, sub_df$feature_id)
message(sprintf("subsampled %d features for %d permutations", length(feats), N_PERMS))

#' Share of the subsample reaching padj < 0.05 under one labelling.
q_hits <- function(labels) {
  p <- vapply(by_feat, function(s) {
    s$model_f <- droplevels(factor(labels[s$study_id], levels = SUBGROUPS))
    if (nlevels(s$model_f) < 2) return(NA_real_)
    fit <- tryCatch(rma(yi = effect_size, sei = se, mods = ~ model_f,
                        data = s, method = METHOD),
                    error = function(e) NULL)
    if (is.null(fit)) NA_real_ else fit$QMp
  }, numeric(1))
  padj <- p.adjust(p, method = "BH")
  c(tested = sum(!is.na(p)), hits = sum(padj < 0.05, na.rm = TRUE))
}

obs <- q_hits(unit_model)
message(sprintf("observed: %d/%d = %.2f%%", obs["hits"], obs["tested"],
                100 * obs["hits"] / obs["tested"]))

#' Each permuted row records its seed, so runs split across processes can be
#' pooled (scripts/pool_permutation_runs.py) and checked for a repeated shuffle.
rows <- data.frame(perm = 0L, seed = NA_integer_, tested = obs["tested"],
                   hits = obs["hits"], pct = 100 * obs["hits"] / obs["tested"],
                   row.names = NULL)
for (p in seq_len(N_PERMS)) {
  set.seed(SEED * 1000L + p)
  perm <- unit_model
  # shuffle labels within species, preserving every stratum's size
  for (sp in unique(unit_species)) {
    idx <- which(unit_species == sp)
    perm[idx] <- sample(unit_model[idx])
  }
  r <- q_hits(perm)
  rows <- rbind(rows, data.frame(perm = p, seed = SEED * 1000L + p,
                                 tested = r["tested"], hits = r["hits"],
                                 pct = 100 * r["hits"] / r["tested"],
                                 row.names = NULL))
  message(sprintf("  perm %2d: %d/%d = %.2f%%", p, r["hits"], r["tested"],
                  100 * r["hits"] / r["tested"]))
}

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write.csv(rows, OUT, row.names = FALSE)

null_pct <- rows$pct[rows$perm > 0]
message(sprintf(
  "\nobserved %.2f%% vs permuted null %.2f%% (range %.2f-%.2f, n=%d)",
  rows$pct[1], mean(null_pct), min(null_pct), max(null_pct), length(null_pct)))
message(sprintf("empirical p = %.4f",
                (sum(null_pct >= rows$pct[1]) + 1) / (length(null_pct) + 1)))
message(sprintf("-> %s", OUT))
