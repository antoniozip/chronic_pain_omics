#!/usr/bin/env Rscript
#' Redraw the manuscript's forest figure from the tracked results.
#'
#' The figure is produced by `pipeline/06_meta_analysis.R`, but re-running that
#' step for transcriptomics takes about two hours and recomputes nothing this
#' figure needs: the per-study effects and the pooled estimates are already on
#' disk. This draws the panel from them, through the same
#' `write_forest_pdf()` the pipeline calls, so the two cannot diverge.
#'
#' Writes manuscript/figures/{modality}_forest.pdf: one page, the top pooled
#' feature, on a page sized to that panel, and {modality}_funnel.pdf for the
#' same feature. The two are a pair and must show the same gene -- the funnel
#' used to take the most-measured feature instead, so the pair described two
#' different genes and neither filename nor title said so. The supplementary
#' set of the top N features stays where step 06 puts it and is not touched.
#'
#' Usage:
#'     Rscript scripts/rebuild_forest_figure.R
#'     Rscript scripts/rebuild_forest_figure.R --modality proteomics

suppressPackageStartupMessages({
  library(metafor)
  library(yaml)
})

REPO_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))
source(file.path(REPO_ROOT, "R", "forest_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
modality <- if (length(args) >= 2 && args[1] == "--modality") args[2] else "transcriptomics"

cfg <- yaml::read_yaml(file.path(REPO_ROOT, "conf", "analysis", "default.yaml"))
meta_dir <- file.path(REPO_ROOT, "results", "meta", modality)
per_study_csv <- file.path(REPO_ROOT, "results", "per_study", modality,
                           "effects_normalized_all.csv")
figures_dir <- file.path(REPO_ROOT, "manuscript", "figures")

for (path in c(file.path(meta_dir, "pooled_effects.csv"), per_study_csv)) {
  if (!file.exists(path)) {
    stop("missing input: ", path, "\nRun the pipeline for this modality first.")
  }
}

pooled <- read.csv(file.path(meta_dir, "pooled_effects.csv"),
                   stringsAsFactors = FALSE)
top_feature <- top_pooled_features(pooled, 1)
if (length(top_feature) == 0) stop("no pooled feature has an adjusted p-value")
message(sprintf("Top pooled feature: %s", top_feature))

# The per-study table is ~285 MB, and only one feature's rows are wanted. Read
# the columns the panel needs and filter, rather than loading it whole.
effects <- read.csv(per_study_csv, stringsAsFactors = FALSE,
                    colClasses = "character")
keep <- effects$feature_id == top_feature
effects <- effects[keep, c("feature_id", "study_id", "effect_size", "se",
                           "effect_unit"), drop = FALSE]
effects$effect_size <- as.numeric(effects$effect_size)
effects$se <- as.numeric(effects$se)
message(sprintf("%d per-study rows for it", nrow(effects)))

out <- write_forest_pdf(
  file.path(figures_dir, sprintf("%s_forest.pdf", modality)),
  effects, top_feature, method = cfg$meta$method,
  slab_max_chars = cfg$meta$slab_max_chars
)
message(sprintf("Wrote %s", out))

# The funnel of the same feature, drawn here rather than by re-running step 06
# for the same reason the forest is: nothing it needs is recomputed by that run.
funnel_path <- file.path(figures_dir, sprintf("%s_funnel.pdf", modality))
if (nrow(effects) >= 3) {
  pdf(funnel_path, width = 7, height = 6)
  tryCatch({
    res <- metafor::rma(yi = effect_size, sei = se, data = effects,
                        method = cfg$meta$method)
    metafor::funnel(res, main = sprintf("%s: %s", modality, top_feature))
  }, error = function(e) message("  funnel failed: ", e$message))
  dev.off()
  message(sprintf("Wrote %s", funnel_path))
} else {
  message("Fewer than 3 studies for it; no funnel drawn")
}
