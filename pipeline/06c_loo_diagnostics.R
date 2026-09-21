#!/usr/bin/env Rscript
#' Leave-one-out sensitivity analysis and influence diagnostics (Step 2.3).
#'
#' For the top-N features in the main meta-analysis, runs:
#'   1. Leave-one-out (LOO): iterates over studies, removes each, re-pools.
#'      Flags features where removing a single study reverses the sign or
#'      makes the estimate non-significant.
#'   2. Weight analysis: flags features where one study contributes >50%
#'      of the inverse-variance weight.
#'   3. Influence diagnostics via metafor::influence().
#'
#' Outputs:
#'   results/meta/transcriptomics/loo/
#'     loo_top_features.csv       — LOO estimates for top-N features × studies
#'     loo_robustness_flags.csv   — per-feature robustness summary + flags
#'     influence_top50.csv        — Cook's D and hat values for top 50
#'     loo_forest_grid.pdf        — forest of LOO estimates for top 12 features

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(ggplot2)
})

REPO_ROOT   <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))), ".."))
EFFECTS_CSV <- file.path(REPO_ROOT, "results", "per_study", "transcriptomics",
                          "effects_normalized_all.csv")
POOLED_CSV  <- file.path(REPO_ROOT, "results", "meta", "transcriptomics", "pooled_effects.csv")
OUT_DIR     <- file.path(REPO_ROOT, "results", "meta", "transcriptomics", "loo")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METHOD      <- "REML"
TOP_N       <- 100    # features to analyse
LOO_TOP_PDF <- 12     # features to plot in grid

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

per_study <- read.csv(EFFECTS_CSV, stringsAsFactors = FALSE)
per_study <- per_study[!is.na(per_study$se) & per_study$se > 0, ]

pooled <- read.csv(POOLED_CSV, stringsAsFactors = FALSE)
pooled <- pooled[order(pooled$pval_pooled), ]
# Top-N by overall significance, requiring k >= 3 (LOO requires k >= 2)
pooled <- pooled[!is.na(pooled$pval_pooled) & !is.na(pooled$k) & pooled$k >= 3, ]
top_features <- head(pooled$feature_id, TOP_N)

message(sprintf("Running LOO for top %d features (k>=3)", length(top_features)))

# ---------------------------------------------------------------------------
# Leave-one-out function
# ---------------------------------------------------------------------------

run_loo <- function(feat, studies_df) {
  studies <- unique(studies_df$study_id)
  rows <- list()
  for (excl in studies) {
    sub <- studies_df[studies_df$study_id != excl, ]
    sub <- sub[!is.na(sub$se) & sub$se > 0, ]
    if (nrow(sub) < 2) next
    res <- tryCatch(
      rma(yi = effect_size, sei = se, data = sub, method = METHOD),
      error = function(e) NULL
    )
    if (is.null(res)) next
    rows[[length(rows) + 1]] <- data.frame(
      feature_id   = feat,
      excl_study   = excl,
      k_remaining  = nrow(sub),
      yi_loo       = res$b[1],
      se_loo       = res$se,
      pval_loo     = res$pval,
      ci_lb_loo    = res$ci.lb,
      ci_ub_loo    = res$ci.ub,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Run LOO for all top features
# ---------------------------------------------------------------------------

loo_all <- list()

for (feat in top_features) {
  sub <- per_study[per_study$feature_id == feat, ]
  loo_res <- run_loo(feat, sub)
  if (!is.null(loo_res) && nrow(loo_res) > 0) {
    loo_all[[feat]] <- loo_res
  }
}

loo_df <- do.call(rbind, Filter(Negate(is.null), loo_all))
write.csv(loo_df, file.path(OUT_DIR, "loo_top_features.csv"), row.names = FALSE)
message(sprintf("LOO table: %d rows", nrow(loo_df)))

# ---------------------------------------------------------------------------
# Robustness flags
# ---------------------------------------------------------------------------

flags <- pooled[pooled$feature_id %in% top_features, ] %>%
  select(feature_id, k, pval_pooled, yi_pooled, padj_pooled) %>%
  mutate(
    n_sign_reversals = NA_integer_,  # studies whose removal reverses sign
    n_sig_losses     = NA_integer_,  # studies whose removal loses significance
    max_single_weight_pct = NA_real_,  # max weight % of one study
    robust           = NA
  )

for (i in seq_len(nrow(flags))) {
  feat     <- flags$feature_id[i]
  main_yi  <- flags$yi_pooled[i]
  main_sig <- !is.na(flags$pval_pooled[i]) && flags$pval_pooled[i] < 0.05

  loo_sub  <- loo_df[loo_df$feature_id == feat, ]
  if (nrow(loo_sub) == 0) next

  flags$n_sign_reversals[i] <- sum(sign(loo_sub$yi_loo) != sign(main_yi), na.rm = TRUE)
  if (main_sig) {
    flags$n_sig_losses[i] <- sum(loo_sub$pval_loo >= 0.05, na.rm = TRUE)
  } else {
    flags$n_sig_losses[i] <- 0L
  }

  # Weight concentration
  feat_sub <- per_study[per_study$feature_id == feat & !is.na(per_study$se) & per_study$se > 0, ]
  if (nrow(feat_sub) > 1) {
    w <- 1 / feat_sub$se^2
    flags$max_single_weight_pct[i] <- max(w) / sum(w) * 100
  }

  flags$robust[i] <- flags$n_sign_reversals[i] == 0 &&
    (is.na(flags$n_sig_losses[i]) || flags$n_sig_losses[i] == 0) &&
    (is.na(flags$max_single_weight_pct[i]) || flags$max_single_weight_pct[i] < 50)
}

flags <- flags %>% arrange(pval_pooled)
write.csv(flags, file.path(OUT_DIR, "loo_robustness_flags.csv"), row.names = FALSE)

robust_n    <- sum(flags$robust == TRUE, na.rm = TRUE)
fragile_n   <- sum(flags$n_sign_reversals > 0, na.rm = TRUE)
weight_conc <- sum(flags$max_single_weight_pct >= 50, na.rm = TRUE)

message(sprintf(
  "\nRobustness summary:\n  Robust (no reversals, no sig loss, weight <50%%): %d/%d\n  Sign reversal in any LOO: %d\n  Weight-concentrated (>50%% in 1 study): %d",
  robust_n, nrow(flags), fragile_n, weight_conc))

# ---------------------------------------------------------------------------
# Influence diagnostics for top-50
# ---------------------------------------------------------------------------

inf_rows <- list()
top50 <- head(flags$feature_id[!is.na(flags$pval_pooled)], 50)

for (feat in top50) {
  sub <- per_study[per_study$feature_id == feat & !is.na(per_study$se) & per_study$se > 0, ]
  if (nrow(sub) < 3) next
  res <- tryCatch(
    rma(yi = effect_size, sei = se, data = sub, method = METHOD, slab = sub$study_id),
    error = function(e) NULL
  )
  if (is.null(res)) next
  inf <- tryCatch(influence(res), error = function(e) NULL)
  if (is.null(inf)) next
  inf_df <- as.data.frame(inf$inf)
  inf_df$feature_id <- feat
  inf_df$study_id   <- sub$study_id
  inf_rows[[feat]] <- inf_df
}

if (length(inf_rows) > 0) {
  inf_combined <- do.call(rbind, inf_rows)
  write.csv(inf_combined, file.path(OUT_DIR, "influence_top50.csv"), row.names = FALSE)
  message(sprintf("Influence table: %d rows", nrow(inf_combined)))
}

# ---------------------------------------------------------------------------
# LOO forest grid for top-12 features
# ---------------------------------------------------------------------------

top12 <- head(flags$feature_id, LOO_TOP_PDF)
pdf(file.path(OUT_DIR, "loo_forest_grid.pdf"), width = 14, height = 3.5 * ceiling(LOO_TOP_PDF / 3))
par(mfrow = c(ceiling(LOO_TOP_PDF / 3), 3), mar = c(3, 4, 2, 1))

for (feat in top12) {
  loo_sub <- loo_df[loo_df$feature_id == feat, ]
  if (nrow(loo_sub) == 0) next

  main_row <- flags[flags$feature_id == feat, ]
  main_yi  <- main_row$yi_pooled
  main_k   <- main_row$k

  xlim_val <- range(c(loo_sub$ci_lb_loo, loo_sub$ci_ub_loo,
                      main_yi - 0.5, main_yi + 0.5), na.rm = TRUE)

  # Color: red if reversal
  cols <- ifelse(sign(loo_sub$yi_loo) != sign(main_yi), "#d6604d", "#2166ac")

  plot(loo_sub$yi_loo, seq_len(nrow(loo_sub)),
       xlim = xlim_val,
       ylim = c(0, nrow(loo_sub) + 2),
       xlab = "log2FC (LOO pooled)",
       yaxt = "n", ylab = "",
       pch  = 16, col = cols,
       main = sprintf("%s  (k=%d)", feat, main_k),
       cex.main = 0.85, cex = 0.8)

  arrows(loo_sub$ci_lb_loo, seq_len(nrow(loo_sub)),
         loo_sub$ci_ub_loo, seq_len(nrow(loo_sub)),
         angle = 90, code = 3, length = 0.03, col = cols, lwd = 0.6)
  axis(2, at = seq_len(nrow(loo_sub)), labels = loo_sub$excl_study,
       las = 2, cex.axis = 0.55)
  abline(v = 0, lty = 2, col = "grey")
  abline(v = main_yi, lty = 1, col = "black", lwd = 1.2)
}

dev.off()
message(sprintf("LOO forest grid → %s", file.path(OUT_DIR, "loo_forest_grid.pdf")))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat("\n========== LOO Robustness Summary ==========\n")
cat(sprintf("Top %d features analysed (k>=3)\n", nrow(flags)))
cat(sprintf("  Robust (no reversals, no sig loss, weight <50%%): %d\n", robust_n))
cat(sprintf("  Features with sign reversal in any LOO:           %d\n", fragile_n))
cat(sprintf("  Features with >50%% weight in single study:        %d\n", weight_conc))
cat("\nTop 10 robust features:\n")
robust_top <- flags[flags$robust == TRUE, ][1:min(10, sum(flags$robust == TRUE, na.rm=TRUE)), ]
print(robust_top[, c("feature_id","k","pval_pooled","yi_pooled","max_single_weight_pct")])
cat("\nTop 5 fragile (sign-reversal) features:\n")
fragile <- flags[!is.na(flags$n_sign_reversals) & flags$n_sign_reversals > 0, ]
fragile <- fragile[order(-fragile$n_sign_reversals), ]
if (nrow(fragile) > 0) {
  print(head(fragile[, c("feature_id","k","pval_pooled","n_sign_reversals","max_single_weight_pct")], 5))
}

message("\nStep 06c complete. Output: ", OUT_DIR)
