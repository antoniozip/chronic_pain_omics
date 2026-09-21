#!/usr/bin/env Rscript
#' Random-effects meta-analysis pipeline step.
#'
#' Reads all per-study effect-size CSVs from results/per_study/{modality}/ and
#' runs a random-effects meta-analysis using metafor (REML by default).
#'
#' Output per modality:
#'   results/meta/{modality}/pooled_effects.csv  — pooled estimates
#'   results/meta/{modality}/heterogeneity.csv   — I², τ², Q, Cochran's Q p-value
#'   manuscript/figures/{modality}_forest.pdf    — forest plot (top N features)
#'   manuscript/figures/{modality}_funnel.pdf    — funnel plot (for bias assessment)
#'
#' Usage:
#'   Rscript pipeline/06_meta_analysis.R --modality transcriptomics
#'   Rscript pipeline/06_meta_analysis.R --modality all --method REML

suppressPackageStartupMessages({
  library(optparse)
  library(metafor)
  library(yaml)
  library(dplyr)
  library(ggplot2)
})

script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
REPO_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
ALL_MODALITIES <- c("genomics", "transcriptomics", "proteomics", "metabolomics", "lipidomics")

source(file.path(REPO_ROOT, "R", "meta_utils.R"))
source(file.path(REPO_ROOT, "R", "forest_utils.R"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

option_list <- list(
  make_option("--modality",        type = "character", default = "transcriptomics"),
  make_option("--analysis-config", type = "character",
              default = file.path(REPO_ROOT, "conf", "analysis", "default.yaml")),
  make_option("--per-study-dir",   type = "character",
              default = file.path(REPO_ROOT, "results", "per_study")),
  make_option("--meta-dir",        type = "character",
              default = file.path(REPO_ROOT, "results", "meta")),
  make_option("--figures-dir",     type = "character",
              default = file.path(REPO_ROOT, "manuscript", "figures")),
  make_option("--method",          type = "character", default = NULL,
              help = "metafor method: REML|DL|ML|HE (overrides config)"),
  make_option("--dry-run",         action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))
cfg <- read_yaml(opt[["analysis-config"]])

`%||%` <- function(a, b) if (!is.null(a)) a else b

method <- opt$method %||% cfg$meta$method
modalities <- if (opt$modality == "all") ALL_MODALITIES else opt$modality

# `single_cell` is poolable but is NOT in ALL_MODALITIES, so `--modality all`
# cannot reach it. The pseudobulked single-cell arm is a stratum of its own and
# is never merged into the bulk transcriptomic pool; naming it is the only way
# to analyse it, and the guard below turns a typo into an error rather than an
# empty output directory.
EXTRA_MODALITIES <- c("single_cell")
unknown <- setdiff(modalities, c(ALL_MODALITIES, EXTRA_MODALITIES))
if (length(unknown) > 0) {
  stop(sprintf("unknown modality: %s", paste(unknown, collapse = ", ")))
}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

load_per_study_effects <- function(per_study_dir, modality) {
  mod_dir <- file.path(per_study_dir, modality)
  # Three tiers, most-processed first, one file per accession:
  #   _effects_gene.csv       re-keyed onto gene symbols from another feature
  #                           space (05i, peptides -> genes for PXD013362)
  #   _effects_normalized.csv gene-symbol-level (05b)
  #   _effects.csv            raw, native identifiers
  # The gene tier exists because a study can be gene-keyed without having gone
  # through 05b's probe normalization. For PXD013362 the _normalized file is
  # peptide-level and enters the pool as 14 pseudo-studies, so preferring the
  # gene file is also what keeps that study's samples independent; see
  # plan/2026-08-14-proteomics-expansion.md.
  # These globs are recursive, and `results/per_study/transcriptomics` holds a
  # pre-canonicalisation backup of 14 studies. It is skipped only because its
  # name starts with a dot and `list.files` defaults to `all.files = FALSE`.
  # Renaming that directory without a dot would silently load those 14 studies
  # a second time, with their stale uncanonicalised symbols.
  gene_files <- list.files(mod_dir, pattern = "_effects_gene\\.csv$", full.names = TRUE,
                           recursive = TRUE)
  norm_files <- list.files(mod_dir, pattern = "_effects_normalized\\.csv$", full.names = TRUE,
                           recursive = TRUE)
  raw_files  <- list.files(mod_dir, pattern = "_effects\\.csv$", full.names = TRUE, recursive = TRUE)
  raw_files  <- raw_files[!grepl("_normalized", raw_files)]

  gene_acc <- sub("_effects_gene\\.csv$", "", basename(gene_files))
  norm_acc <- sub("_effects_normalized\\.csv$", "", basename(norm_files))
  raw_acc  <- sub("_effects\\.csv$", "", basename(raw_files))
  norm_files <- norm_files[!norm_acc %in% gene_acc]
  norm_acc   <- norm_acc[!norm_acc %in% gene_acc]
  only_raw <- raw_files[!raw_acc %in% c(gene_acc, norm_acc)]

  files_to_load <- c(gene_files, norm_files, only_raw)
  # Drop studies superseded by a finer-grained re-analysis of the same samples.
  # This loader globs the per-study tables directly, so filtering the combined
  # effects_normalized_all.csv does not reach it.
  superseded <- superseded_studies(
    file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv"))
  if (length(superseded) > 0) {
    accs <- sub("_effects(_normalized|_gene)?\\.csv$", "", basename(files_to_load))
    drop <- accs %in% superseded
    for (acc in accs[drop]) {
      message(sprintf("[%s] excluding %s: superseded (conf/analysis/superseded_studies.csv)",
                      modality, acc))
    }
    files_to_load <- files_to_load[!drop]
    gene_files <- gene_files[!sub("_effects_gene\\.csv$", "", basename(gene_files)) %in% superseded]
    norm_files <- norm_files[!sub("_effects_normalized\\.csv$", "", basename(norm_files)) %in% superseded]
    only_raw   <- only_raw[!sub("_effects\\.csv$", "", basename(only_raw)) %in% superseded]
  }
  if (length(files_to_load) == 0) {
    message(sprintf("[%s] No effect files found in %s", modality, mod_dir))
    return(NULL)
  }
  df <- bind_rows(lapply(files_to_load, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$feature_id <- as.character(d$feature_id)
    if ("original_feature_id" %in% colnames(d))
      d$original_feature_id <- as.character(d$original_feature_id)
    d
  }))
  message(sprintf("[%s] Loaded %d effect rows from %d studies (%d gene-keyed, %d normalized, %d raw)",
                  modality, nrow(df), length(files_to_load),
                  length(gene_files), length(norm_files), length(only_raw)))
  df
}


save_results <- function(pooled, hetero, modality, meta_dir) {
  dir.create(file.path(meta_dir, modality), recursive = TRUE, showWarnings = FALSE)
  write.csv(pooled, file.path(meta_dir, modality, "pooled_effects.csv"), row.names = FALSE)
  write.csv(hetero, file.path(meta_dir, modality, "heterogeneity.csv"), row.names = FALSE)
  message(sprintf("[%s] Results written to %s/%s/", modality, meta_dir, modality))
}


# ---------------------------------------------------------------------------
# Meta-analysis per feature
# ---------------------------------------------------------------------------

meta_one_feature <- function(feature_df, method) {
  feature_df <- feature_df[!is.na(feature_df$se) & feature_df$se > 0, ]
  if (nrow(feature_df) == 0) return(NULL)

  study_ids <- paste(unique(feature_df$study_id), collapse = ";")

  # Single-study pass-through: record as k=1, no pooling
  if (nrow(feature_df) == 1) {
    row <- feature_df[1, ]
    return(data.frame(
      feature_id    = row$feature_id,
      k             = 1L,
      study_ids     = study_ids,
      yi_pooled     = row$effect_size,
      se_pooled     = row$se,
      ci_lb         = row$effect_size - 1.96 * row$se,
      ci_ub         = row$effect_size + 1.96 * row$se,
      pval_pooled   = row$pval,
      I2            = NA_real_,
      tau2          = NA_real_,
      Q             = NA_real_,
      Q_pval        = NA_real_,
      H2            = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  res <- tryCatch(
    rma(yi = effect_size, sei = se, data = feature_df,
        method = method, slab = feature_df$study_id),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  data.frame(
    feature_id    = feature_df$feature_id[1],
    k             = res$k,
    study_ids     = study_ids,
    yi_pooled     = res$b[1],
    se_pooled     = res$se,
    ci_lb         = res$ci.lb,
    ci_ub         = res$ci.ub,
    pval_pooled   = res$pval,
    I2            = res$I2,
    tau2          = res$tau2,
    Q             = res$QE,
    Q_pval        = res$QEp,
    H2            = res$H2,
    stringsAsFactors = FALSE
  )
}


run_meta_analysis <- function(effects_df, method, cfg) {
  features <- split(effects_df, effects_df$feature_id)
  message(sprintf("  Running meta-analysis for %d features across %d studies ...",
                  length(features),
                  length(unique(effects_df$study_id))))

  results <- lapply(features, meta_one_feature, method = method)
  pooled  <- bind_rows(Filter(Negate(is.null), results))
  if (nrow(pooled) == 0) return(list(pooled = pooled, heterogeneity = data.frame()))

  pooled <- pooled %>%
    arrange(pval_pooled) %>%
    mutate(heterogeneity_level = case_when(
             I2 >= cfg$meta$i2_high     ~ "high",
             I2 >= cfg$meta$i2_moderate ~ "moderate",
             TRUE                       ~ "low"
           ))
  # k=1 features are pass-throughs of a single study's raw p-value, not
  # meta-analytic ones; they stay in the output with padj_pooled = NA and are
  # excluded from the correction family. See R/meta_utils.R.
  pooled <- adjust_pooled_fdr(pooled, pval_col = "pval_pooled",
                              min_k = cfg$meta$min_studies %||% 2L,
                              out_col = "padj_pooled")

  hetero_summary <- pooled %>%
    summarise(
      n_features    = n(),
      n_pooled      = sum(k >= (cfg$meta$min_studies %||% 2L), na.rm = TRUE),
      n_sig         = sum(padj_pooled < 0.05, na.rm = TRUE),
      median_I2     = median(I2, na.rm = TRUE),
      median_tau2   = median(tau2, na.rm = TRUE),
      pct_high_het  = 100 * mean(heterogeneity_level == "high", na.rm = TRUE)
    )

  list(pooled = pooled, heterogeneity = hetero_summary)
}


# ---------------------------------------------------------------------------
# Publication bias (Egger's test)
# ---------------------------------------------------------------------------

run_egger_test <- function(effects_df, modality) {
  features_2plus <- effects_df %>%
    filter(!is.na(se), se > 0) %>%
    group_by(feature_id) %>%
    filter(n() >= 10) %>%        # Egger's requires ≥10 studies
    ungroup()

  if (nrow(features_2plus) == 0) {
    message(sprintf("  [%s] Too few multi-study features for Egger's test", modality))
    return(invisible(NULL))
  }

  # rma() resolves yi/sei by non-standard evaluation and cannot see dplyr's
  # data mask, so calling it inside summarise() raised "Cannot find the
  # object/variable ('effect_size')" for every feature. The errors were
  # swallowed into NA and filtered out, so the test reported 0 / 0 while
  # never evaluating anything. Each feature is now fitted once, with data=
  # supplied explicitly, and the single fit serves both statistics.
  per_feature <- split(features_2plus, features_2plus$feature_id)
  egger_results <- bind_rows(lapply(names(per_feature), function(feat) {
    sub <- per_feature[[feat]]
    res <- tryCatch({
      fit <- rma(yi = effect_size, sei = se, data = sub, method = "FE")
      reg <- regtest(fit)
      data.frame(feature_id = feat, k = nrow(sub),
                 egger_z = reg$zval, egger_p = reg$pval,
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    res
  }))

  n_failed <- length(per_feature) - nrow(egger_results)
  if (n_failed > 0) {
    message(sprintf("  [%s] Egger's test: %d of %d features could not be fitted",
                    modality, n_failed, length(per_feature)))
  }
  if (nrow(egger_results) == 0) {
    message(sprintf("  [%s] Egger's test produced no results", modality))
    return(invisible(NULL))
  }

  egger_results <- egger_results %>%
    mutate(bias_flagged = egger_p < cfg$meta$egger_alpha)

  n_flagged <- sum(egger_results$bias_flagged, na.rm = TRUE)
  message(sprintf("  [%s] Egger's test: %d / %d features flagged for potential bias",
                  modality, n_flagged, nrow(egger_results)))
  egger_results
}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

plot_forest <- function(effects_df, pooled_df, modality, figures_dir, cfg) {
  top_features <- top_pooled_features(pooled_df, cfg$meta$top_n_features)
  if (length(top_features) == 0) return(invisible(NULL))

  # Two files, because a pdf() device fixes one page size for the whole
  # document and the two readers need different ones. The supplementary set
  # is sized for its deepest panel, which is right for a file read on screen;
  # the manuscript figure shows only the top feature, so it gets a page the
  # size of that panel. Drawing the manuscript figure as page 1 of the
  # supplementary set gave it a 10x30 inch page for a 3-study forest, which
  # LaTeX refused to place and which blocked every float behind it.
  supp_path <- file.path(figures_dir,
                         sprintf("%s_forest_supplementary.pdf", modality))
  write_forest_pdf(supp_path, effects_df, top_features, method = method,
                   slab_max_chars = cfg$meta$slab_max_chars)
  message(sprintf("  Forest plots (supplementary) → %s", supp_path))

  main_path <- file.path(figures_dir, sprintf("%s_forest.pdf", modality))
  write_forest_pdf(main_path, effects_df, top_features[1], method = method,
                   slab_max_chars = cfg$meta$slab_max_chars)
  message(sprintf("  Forest plot (manuscript) → %s", main_path))
}


plot_funnel <- function(effects_df, pooled_df, modality, figures_dir) {
  # The same feature the manuscript forest draws, which is the first of
  # top_pooled_features() -- ranked within the BH-adjusted family, never across
  # all rows. The funnel used to take the most-measured feature instead, so the
  # two figures of a pair described different genes and neither filename nor
  # caption said so: the single-cell arm's forest showed Ppm1l while its funnel
  # showed C3. Small-study behaviour is only worth reading for the feature the
  # paper is actually making a claim about.
  top_feat <- top_pooled_features(pooled_df, 1)
  if (length(top_feat) == 0) return(invisible(NULL))
  top_feat <- top_feat[[1]]
  feat_df  <- effects_df[effects_df$feature_id == top_feat & !is.na(effects_df$se), ]
  if (nrow(feat_df) < 3) return(invisible(NULL))

  out_path <- file.path(figures_dir, sprintf("%s_funnel.pdf", modality))
  pdf(out_path, width = 7, height = 6)
  tryCatch({
    res <- rma(yi = effect_size, sei = se, data = feat_df, method = method)
    # ASCII separator: the pdf() device cannot encode an em-dash in its
    # default font and silently substitutes dots, so the title read
    # "single_cell ... Ppm1l".
    funnel(res, main = sprintf("%s: %s", modality, top_feat))
  }, error = function(e) NULL)
  dev.off()
  message(sprintf("  Funnel plot → %s", out_path))
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

for (modality in modalities) {
  message("\n========== ", toupper(modality), " ==========")

  effects_df <- load_per_study_effects(opt[["per-study-dir"]], modality)
  if (is.null(effects_df) || nrow(effects_df) == 0) next

  if (opt[["dry-run"]]) {
    message(sprintf("[dry-run][%s] Would meta-analyse %d features from %d studies",
                    modality,
                    length(unique(effects_df$feature_id)),
                    length(unique(effects_df$study_id))))
    next
  }

  ma_result <- run_meta_analysis(effects_df, method, cfg)
  pooled  <- ma_result$pooled
  hetero  <- ma_result$heterogeneity

  if (nrow(pooled) == 0) {
    message(sprintf("[%s] No poolable features.", modality))
    next
  }

  save_results(pooled, hetero, modality, opt[["meta-dir"]])
  # The return value was previously discarded, so even a working test left no
  # artifact behind. Publication-bias results are now written beside the
  # pooled estimates.
  egger <- run_egger_test(effects_df, modality)
  if (!is.null(egger) && nrow(egger) > 0) {
    egger_path <- file.path(opt[["meta-dir"]], modality, "egger_test.csv")
    write.csv(egger, egger_path, row.names = FALSE)
    message(sprintf("[%s] Egger results -> %s", modality, egger_path))
  }

  dir.create(opt[["figures-dir"]], recursive = TRUE, showWarnings = FALSE)
  plot_forest(effects_df, pooled, modality, opt[["figures-dir"]], cfg)
  plot_funnel(effects_df, pooled, modality, opt[["figures-dir"]])

  # Rank within the adjusted (k>=2) family so single-study pass-throughs are
  # never presented as the modality's leading meta-analytic signals.
  ranked <- pooled[!is.na(pooled$padj_pooled), ]
  message(sprintf("[%s] Top 5 pooled features (k>=2) by adjusted p-value:", modality))
  print(head(ranked[, c("feature_id", "k", "yi_pooled", "I2", "padj_pooled")], 5))
}

message("\nStep 06 complete.")
