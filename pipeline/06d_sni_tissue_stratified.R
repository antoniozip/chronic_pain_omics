#!/usr/bin/env Rscript
#' SNI Tissue Analysis — Fast, focused on top features.
#' Usage: Rscript --no-init-file pipeline/06d_sni_tissue_stratified.R

suppressPackageStartupMessages({
  library(metafor)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))), ".."))
EFFECTS_CSV <- file.path(REPO_ROOT, "results", "per_study", "transcriptomics",
                          "effects_normalized_all.csv")
TISSUE_CSV  <- file.path(REPO_ROOT, "data", "interim", "sni_tissue_annotations.csv")
OUT_DIR     <- file.path(REPO_ROOT, "results", "meta", "transcriptomics", "stratified")

# moderator_df() lives in R/meta_utils.R. metafor's `fit$m` is the number of
# moderator coefficients and `fit$p` is that count plus the intercept; neither
# is the quantity to publish or to gate on, `QMdf[1]` is.
source(file.path(REPO_ROOT, "R", "meta_utils.R"))


#' The pre-specified tissue-moderator family.
#'
#' Every feature observed in at least `min_tissues` tissues across at least
#' `min_studies` studies. The rule is a property of the stratum's design and is
#' evaluated before any effect or p-value is read, which is what makes the BH
#' adjustment over this family valid.
#'
#' This used to be `head(x, 200)` after `order(-abs(yi))` — the family was
#' selected on effect magnitude, so the retained 200 were enriched for small
#' p-values and the adjusted count was anti-conservative with respect to
#' correcting over all candidate features. A cap is no longer needed: the
#' features are now indexed once by `split()` instead of re-scanning the whole
#' table per feature.
#'
#' @param df_st Single-tissue rows of the stratum.
#' @return Character vector of feature ids (deterministic, alphabetical).
sni_tissue_family <- function(df_st, min_tissues = 2L, min_studies = 3L) {
  n_tissues <- tapply(df_st$tissue, df_st$feature_id, function(x) length(unique(x)))
  n_studies <- tapply(df_st$study_id, df_st$feature_id, function(x) length(unique(x)))
  ids <- names(n_tissues)
  ids[!is.na(n_tissues[ids]) & !is.na(n_studies[ids]) &
        n_tissues[ids] >= min_tissues & n_studies[ids] >= min_studies]
}


#' Fit the tissue moderator for one feature.
#'
#' @param sub Rows of one feature (single-tissue studies only).
#' @param fid Feature id, written into the result row.
#' @return One-row data.frame, or NULL when the test is not identifiable.
sni_tissue_moderator <- function(sub, fid) {
  sub <- sub[!is.na(sub$se) & sub$se > 0, , drop = FALSE]
  tissues <- unique(sub$tissue)
  if (length(tissues) < 2 || nrow(sub) < 3) return(NULL)

  fit <- tryCatch(
    metafor::rma(yi = effect_size, sei = se, mods = ~ tissue, data = sub,
                 method = "REML"),
    error = function(e) NULL
  )
  # At least one moderator df. The gate used to be `fit$m < 2`, which under
  # metafor's definition of `$m` (moderator coefficients, intercept excluded)
  # silently skipped every two-tissue feature even though the family rule
  # admits them.
  if (is.null(fit) || moderator_df(fit) < 1) return(NULL)

  data.frame(
    feature_id = fid, k = fit$k, n_tissues = length(tissues),
    I2_residual = fit$I2, tau2_residual = fit$tau2,
    Q_between = fit$QM, Q_between_df = moderator_df(fit),
    Q_between_pval = fit$QMp,
    Q_residual = fit$QE, Q_residual_pval = fit$QEp,
    stringsAsFactors = FALSE
  )
}


# Load tissue map
tissue_map <- read.csv(TISSUE_CSV, stringsAsFactors = FALSE)

# Load per-study effects for SNI only
df <- read.csv(EFFECTS_CSV, stringsAsFactors = FALSE)
df <- df[!is.na(df$se) & df$se > 0, ]
df <- df[df$study_id %in% tissue_map$study_id, ]
df$tissue <- tissue_map$tissue_category[match(df$study_id, tissue_map$study_id)]

# Only single-tissue studies
df_st <- df[df$tissue != "multi_tissue", ]

cat(sprintf("SNI single-tissue: %d rows, %d studies, %d features\n",
            nrow(df_st), length(unique(df_st$study_id)), length(unique(df_st$feature_id))))

# Tissue distribution
cat("Tissues:", paste(sort(unique(df_st$tissue)), collapse=", "), "\n")

# ---------------------------------------------------------------------------
# Tissue moderator test over the pre-specified family
# ---------------------------------------------------------------------------

family <- sni_tissue_family(df_st)
cat(sprintf("Features in ≥2 tissues & ≥3 studies (pre-specified family): %d\n",
            length(family)))

# Index once: a linear scan of the stratum per feature is what the old
# top-200 cap was really buying.
sub_by_feature <- split(df_st, df_st$feature_id)

results <- list()
for (fid in family) {
  row <- sni_tissue_moderator(sub_by_feature[[fid]], fid)
  if (!is.null(row)) results[[fid]] <- row
}

meta_df <- do.call(rbind, results)
if (!is.null(meta_df) && nrow(meta_df) > 0) {
  meta_df$padj <- p.adjust(meta_df$Q_between_pval, method = "BH")
  meta_df <- meta_df[order(meta_df$padj), ]
  
  n_sig <- sum(meta_df$padj < 0.05, na.rm = TRUE)
  cat(sprintf("\n=== Results: %d genes tested, %d significant (padj<0.05) ===\n\n",
              nrow(meta_df), n_sig))
  
  print(meta_df[, c("feature_id","k","n_tissues","Q_between","Q_between_pval","padj","I2_residual")],
        row.names = FALSE)
  
  write.csv(meta_df, file.path(OUT_DIR, "sni_tissue_meta.csv"), row.names = FALSE)
  
  # Per-tissue effects for top hits
  top_hits <- head(meta_df$feature_id, 15)
  per_tissue <- list()
  for (fid in top_hits) {
    sub <- df_st[df_st$feature_id == fid, ]
    for (tis in unique(sub$tissue)) {
      sub_t <- sub[sub$tissue == tis, ]
      fit_t <- tryCatch({
        if (nrow(sub_t) == 1) {
          list(b = sub_t$effect_size[1], se = sub_t$se[1],
               ci.lb = sub_t$effect_size[1] - 1.96*sub_t$se[1],
               ci.ub = sub_t$effect_size[1] + 1.96*sub_t$se[1],
               k = 1)
        } else {
          rma(yi = effect_size, sei = se, data = sub_t, method = "REML")
        }
      }, error = function(e) NULL)
      if (!is.null(fit_t)) {
        per_tissue[[paste(fid, tis, sep="|")]] <- data.frame(
          feature_id = fid, tissue = tis,
          k = fit_t$k, yi = as.numeric(fit_t$b[1]), se = fit_t$se[1],
          ci_lb = fit_t$ci.lb[1], ci_ub = fit_t$ci.ub[1],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pt_df <- do.call(rbind, per_tissue)
  if (!is.null(pt_df)) {
    write.csv(pt_df, file.path(OUT_DIR, "sni_per_tissue_effects.csv"), row.names = FALSE)
    cat(sprintf("Saved %d per-tissue estimates\n", nrow(pt_df)))
  }
} else {
  cat("No genes passed tissue moderator test.\n")
  meta_df <- data.frame()
}

cat("\nDone.\n")
