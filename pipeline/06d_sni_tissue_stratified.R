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

# Load SNI pooled results
sni_pooled <- read.csv(file.path(OUT_DIR, "SNI_pooled.csv"))

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
# Focus on top SNI features (by pooled |yi|) with ≥2 tissues
# ---------------------------------------------------------------------------

# Get top features: those with k>=2 and present in ≥2 tissues
feat_tissues <- split(df_st$tissue, df_st$feature_id) |>
  sapply(function(x) length(unique(x)))
feat_studies <- split(df_st$study_id, df_st$feature_id) |>
  sapply(function(x) length(unique(x)))

multi_tissue <- names(feat_tissues)[feat_tissues >= 2 & feat_studies >= 3]
cat(sprintf("Features in ≥2 tissues & ≥3 studies: %d\n", length(multi_tissue)))

# Merge SNI pooled effect sizes to prioritize
sni_pooled_idx <- sni_pooled[match(multi_tissue, sni_pooled$feature_id), ]
sni_pooled_idx <- sni_pooled_idx[!is.na(sni_pooled_idx$yi), ]
sni_pooled_idx <- sni_pooled_idx[order(-abs(sni_pooled_idx$yi)), ]

# Top 200 features by effect magnitude
top_features <- head(sni_pooled_idx$feature_id, 200)
cat(sprintf("Testing top %d features by |yi|\n", length(top_features)))

# ---------------------------------------------------------------------------
# Fast tissue moderator test (top features only)
# ---------------------------------------------------------------------------

results <- list()
for (i in seq_along(top_features)) {
  fid <- top_features[i]
  sub <- df_st[df_st$feature_id == fid & df_st$tissue != "multi_tissue", ]
  sub <- sub[!is.na(sub$se) & sub$se > 0, ]
  
  tissues <- unique(sub$tissue)
  if (length(tissues) < 2 || nrow(sub) < 3) next
  
  fit <- tryCatch(
    rma(yi = effect_size, sei = se, mods = ~ tissue, data = sub, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$m < 2) next
  
  results[[fid]] <- data.frame(
    feature_id = fid, k = fit$k, n_tissues = length(tissues),
    I2_residual = fit$I2, tau2_residual = fit$tau2,
    Q_between = fit$QM, Q_between_df = fit$m,
    Q_between_pval = fit$QMp,
    Q_residual = fit$QE, Q_residual_pval = fit$QEp,
    stringsAsFactors = FALSE
  )
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
