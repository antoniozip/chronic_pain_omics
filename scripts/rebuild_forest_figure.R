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
#' different genes and neither filename nor title said so.
#'
#' With --supplementary it redraws {modality}_forest_supplementary.pdf, the
#' top N features one per page, instead: the same panels step 06 writes, from
#' the same tables, so the set can be reshaped without the two-hour run.
#'
#' Usage:
#'     Rscript scripts/rebuild_forest_figure.R
#'     Rscript scripts/rebuild_forest_figure.R --modality proteomics
#'     Rscript scripts/rebuild_forest_figure.R --supplementary

suppressPackageStartupMessages({
  library(metafor)
  library(yaml)
})

REPO_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))
source(file.path(REPO_ROOT, "R", "forest_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args); if (!is.na(i) && length(args) > i) args[i + 1] else default
}
modality <- argval("--modality", "transcriptomics")
supplementary <- "--supplementary" %in% args

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
n_features <- if (supplementary) cfg$meta$top_n_features else 1
top_feature <- top_pooled_features(pooled, n_features)
if (length(top_feature) == 0) stop("no pooled feature has an adjusted p-value")
message(sprintf("%d top pooled feature(s): %s", length(top_feature),
                paste(utils::head(top_feature, 3), collapse = ", ")))

# The per-study table is ~285 MB, and only the drawn features' rows are
# wanted. Read the columns the panels need and filter, rather than keeping it.
effects <- read.csv(per_study_csv, stringsAsFactors = FALSE,
                    colClasses = "character")
keep <- effects$feature_id %in% top_feature
effects <- effects[keep, c("feature_id", "study_id", "effect_size", "se",
                           "effect_unit"), drop = FALSE]
effects$effect_size <- as.numeric(effects$effect_size)
effects$se <- as.numeric(effects$se)
message(sprintf("%d per-study rows for it", nrow(effects)))

out <- write_forest_pdf(
  file.path(figures_dir, sprintf("%s_forest%s.pdf", modality,
                                 if (supplementary) "_supplementary" else "")),
  effects, top_feature, method = cfg$meta$method,
  slab_max_chars = cfg$meta$slab_max_chars
)
message(sprintf("Wrote %s", out))

# The funnel belongs to the manuscript figure's single feature; a rebuild of
# the supplementary set neither draws nor invalidates it.
if (supplementary) quit(save = "no")

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
