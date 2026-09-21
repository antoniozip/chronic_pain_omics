#!/usr/bin/env Rscript
#' Collapse PXD013362 across brain regions, with region as a moderator.
#'
#' The study ships one contrast per brain region per condition and its effect
#' table encodes that in `study_id` (PXD013362_PAG_OIH and so on). Feeding those
#' to `06` makes one experiment enter the pool as up to 14 independent studies:
#' the k distribution reached 17 across three studies, and every feature in the
#' k>=3 correction family was within-study replication.
#'
#' Regions within a condition are the same animals, so they are replicates, not
#' studies. Migraine and OIH are separate cohorts, so they remain separate
#' units. This script therefore emits **one row per feature per condition**,
#' pooling the regions with a random-effects model and keeping region as a
#' moderator rather than discarding the information:
#'
#'   effect_size   pooled across regions within the condition
#'   se            of that pooled estimate, so between-region heterogeneity
#'                 widens it rather than being ignored
#'   I2, tau2      between-region heterogeneity
#'   region_qm_p   omnibus test that the effect differs by region
#'   n_regions     regions contributing
#'
#' A feature measured in one region only passes through unchanged, with the
#' moderator columns NA.
#'
#' Output: results/per_study/proteomics/PXD013362/PXD013362_effects_collapsed.csv
#'
#' Usage:
#'   Rscript pipeline/05h_collapse_pxd013362.R

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

REPO_ROOT <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))

IN  <- file.path(REPO_ROOT, "results", "per_study", "proteomics", "PXD013362",
                 "PXD013362_effects.csv")
OUT <- file.path(REPO_ROOT, "results", "per_study", "proteomics", "PXD013362",
                 "PXD013362_effects_collapsed.csv")

if (!file.exists(IN)) stop("no per-study table at ", IN)

df <- read.csv(IN, stringsAsFactors = FALSE)
df <- usable_se(df)
message(sprintf("in: %d rows, %d features, %d study_id values, %d regions, %d conditions",
                nrow(df), length(unique(df$feature_id)), length(unique(df$study_id)),
                length(unique(df$region)), length(unique(df$condition))))

collapse_one <- function(sub) {
  cond <- sub$condition[1]
  base <- data.frame(
    feature_id = sub$feature_id[1],
    study_id   = paste0("PXD013362_", cond),
    condition  = cond,
    region     = "pooled",
    n_regions  = length(unique(sub$region)),
    n_control  = sum(sub$n_control, na.rm = TRUE),
    n_treatment = sum(sub$n_treatment, na.rm = TRUE),
    control_mean = mean(sub$control_mean, na.rm = TRUE),
    treatment_mean = mean(sub$treatment_mean, na.rm = TRUE),
    cohort = "", stringsAsFactors = FALSE
  )
  if (nrow(sub) == 1) {
    return(cbind(base, data.frame(
      effect_size = sub$effect_size[1], se = sub$se[1], pval = sub$pval[1],
      I2 = NA_real_, tau2 = NA_real_, region_qm_p = NA_real_)))
  }

  res <- tryCatch(rma(yi = effect_size, sei = se, data = sub, method = "REML"),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)

  # Region as moderator. droplevels so a feature seen in two regions does not
  # carry every level of the factor, which would leave rank-deficient columns
  # and a degrees-of-freedom count that stops tracking the design
  # (see moderator_df in R/meta_utils.R).
  qm_p <- NA_real_
  if (length(unique(sub$region)) > 1) {
    sub$region_f <- droplevels(factor(sub$region))
    mod <- tryCatch(rma(yi = effect_size, sei = se, mods = ~ region_f,
                        data = sub, method = "REML"),
                    error = function(e) NULL)
    if (!is.null(mod)) qm_p <- mod$QMp
  }

  cbind(base, data.frame(
    effect_size = res$b[1], se = res$se, pval = res$pval,
    I2 = res$I2, tau2 = res$tau2, region_qm_p = qm_p))
}

groups <- split(df, list(df$feature_id, df$condition), drop = TRUE)
out <- bind_rows(lapply(groups, collapse_one))
out$padj <- NA_real_
for (cond in unique(out$condition)) {
  idx <- out$condition == cond
  out$padj[idx] <- p.adjust(out$pval[idx], method = "BH")
}

write.csv(out, OUT, row.names = FALSE)

message(sprintf("out: %d rows, %d features, %d study units (%s)",
                nrow(out), length(unique(out$feature_id)),
                length(unique(out$study_id)),
                paste(unique(out$study_id), collapse = ", ")))
message(sprintf("  region-dependent at qm_p<0.05: %d of %d testable",
                sum(out$region_qm_p < 0.05, na.rm = TRUE),
                sum(!is.na(out$region_qm_p))))
message(sprintf("  significant within condition (padj<0.05): %d",
                sum(out$padj < 0.05, na.rm = TRUE)))
message("\nNote: 06 must read the collapsed table, not the original; the two ",
        "cannot both be present in the glob.")
