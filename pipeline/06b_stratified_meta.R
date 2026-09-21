#!/usr/bin/env Rscript
#' Stratified meta-analysis by pain model.
#'
#' Partitions the transcriptomics studies by pain model (rodent CCI/SNI/SNL/
#' CFA/CIPN and human neuropathic/LBP/nociplastic/in vitro) and runs separate
#' random-effects meta-analyses per subgroup. Also runs a formal Q-test for
#' between-model heterogeneity (meta-regression) on genes present in ≥2 models.
#'
#' BH correction within each stratum covers features seen in at least
#' meta.min_studies studies (conf/analysis/default.yaml); features below that
#' stay in the output with padj = NA. See adjust_pooled_fdr in R/meta_utils.R.
#'
#' Outputs:
#'   results/meta/transcriptomics/stratified/
#'     {model}_pooled.csv        — pooled estimates within each pain model
#'     between_model_Q.csv       — Q_between test per feature
#'     model_heatmap.pdf         — heatmap of top feature effects across models
#'     known_genes_forest.pdf    — per-model forest plots for Jun/Cd44/Atf3/Tspo/Cklf

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(ggplot2)
})

REPO_ROOT    <- normalizePath(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))), ".."))
EFFECTS_CSV  <- file.path(REPO_ROOT, "results", "per_study", "transcriptomics",
                           "effects_normalized_all.csv")
OUT_DIR      <- file.path(REPO_ROOT, "results", "meta", "transcriptomics", "stratified")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(REPO_ROOT, "R", "meta_utils.R"))

METHOD <- "REML"
# Minimum studies for a feature to enter the BH family, read from the shared
# analysis config so 06 and 06b cannot drift apart. At 3, any stratum with
# only two contributing studies (CFA) has an empty family by construction.
MIN_STUDIES <- tryCatch(
  yaml::read_yaml(file.path(REPO_ROOT, "conf", "analysis", "default.yaml"))$meta$min_studies,
  error = function(e) 2L)
if (is.null(MIN_STUDIES)) MIN_STUDIES <- 2L
message(sprintf("BH correction family: features in >= %d studies", MIN_STUDIES))

# ---------------------------------------------------------------------------
# Study → pain model mapping
# ---------------------------------------------------------------------------

STUDY_MODEL <- c(
  GSE172133 = "CCI",  GSE180627 = "CCI",  GSE186237 = "CCI",
  GSE217932 = "CCI",  GSE2636   = "CCI",
  GSE145226 = "SNI",  GSE15041  = "SNI",  GSE18803  = "SNI",
  GSE197233 = "SNI",  GSE214204 = "SNI",  GSE237896 = "SNI",
  GSE89224  = "SNI",  GSE91396  = "SNI",
  GSE10238  = "SNL",  GSE24982  = "SNL",  GSE325961 = "SNL",
  GSE63442  = "SNL",
  GSE147216 = "CFA",  GSE159895 = "CFA",
  GSE30691  = "multi", GSE301439 = "other",
  # Manual annotation batch (2026-05-20)
  GSE113941 = "CIPN", GSE143895 = "CCI",  GSE162284 = "FLIT",
  GSE253183 = "CIPN", GSE289097 = "CFA_orofacial",
  GSE306455 = "SNI",  GSE318938 = "CIPN",
  GSE102937 = "SNI",  GSE241361_DRG = "SNI",  GSE241361_Spinal_cord = "SNI",  GSE256472 = "SNI",
  GSE102721 = "SNI",  GSE117526 = "SNI",  GSE138024 = "SNI",
  GSE160543 = "CIPN",  GSE175760 = "SNI",  GSE224814 = "SNI",
  GSE236754 = "SNI",  GSE245768 = "CCI",  GSE272517 = "other",
  GSE154816 = "CIPN",  GSE283281 = "SNI",  GSE295863 = "other",
  GSE198608 = "SNI",  GSE293520 = "other",  GSE60670 = "SNI",
  GSE265957 = "SNI",
  GSE53860 = "SNL",  GSE53764 = "SNL",
  GSE326825 = "other",  GSE130755 = "other",
  GSE92718 = "CCI",  GSE38038 = "SNL",
  GSE51296 = "other",  GSE51295 = "other",  GSE51294 = "other",
  GSE135080 = "other",
  # Human studies (EI-1 expert decisions 2026-05-20)
  GSE126611 = "neuropathic_human",   # blood; NL-1 (neuropathic pain) vs Control
  GSE177034 = "lbp_human",           # blood; Persistent vs Resolved LBP (t0)
  GSE221921 = "nociplastic_human",   # blood; fibromyalgia vs healthy
  GSE250152 = "neuropathic_human",   # peripheral nerve; Morton's neuroma vs control
  GSE235287 = "invitro_human"        # SH-SY5Y cell line; SMN1-OE vs control plasmid
)

# Strata contributed by the 2026-08-28 search amendment are declared in
# conf/analysis/study_strata.csv rather than in this vector. A hardcoded map
# was tolerable for one search; it makes every retrieval change an edit to
# analysis source, and 46 studies arrived at once. The file is merged over
# STUDY_MODEL below, so a study named in both takes the file's assignment.
STRATA_FILE <- file.path(REPO_ROOT, "conf", "analysis", "study_strata.csv")
extra_strata <- character(0)
if (file.exists(STRATA_FILE)) {
  .sf <- read.csv(STRATA_FILE, stringsAsFactors = FALSE)
  .sf <- .sf[!is.na(.sf$study_id) & nzchar(.sf$study_id), ]
  extra_strata <- stats::setNames(as.character(.sf$stratum), .sf$study_id)
  message(sprintf("study_strata.csv: %d studies across %d strata",
                  length(extra_strata), length(unique(extra_strata))))
}

SUBGROUPS <- unique(c("CCI", "SNI", "SNL", "CFA", "CIPN", "neuropathic_human",
                      "lbp_human", "nociplastic_human", "invitro_human",
                      unname(extra_strata)))

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

df <- read.csv(EFFECTS_CSV, stringsAsFactors = FALSE)
df <- usable_se(df)
# The file wins where both name a study, so a stratum can be corrected without
# editing this script. Merging the two before the lookup keeps that precedence
# while letting a derived unit inherit from its parent accession: STUDY_MODEL
# is keyed on accessions, and GSE180627_S1 is in no map by construction.
STRATUM_MAP <- STUDY_MODEL
if (length(extra_strata) > 0) STRATUM_MAP[names(extra_strata)] <- extra_strata
df$pain_model <- resolve_by_parent(df$study_id, STRATUM_MAP)
.inherited <- unique(df$study_id[!is.na(df$pain_model) &
                                 !(df$study_id %in% names(STRATUM_MAP))])
if (length(.inherited) > 0) {
  message(sprintf("Stratum inherited from parent accession for %d derived unit(s): %s",
                  length(.inherited), paste(sort(.inherited), collapse = ", ")))
}
.unmapped <- unique(df$study_id[is.na(df$pain_model)])
if (length(.unmapped) > 0) {
  message(sprintf("No stratum for %d study(ies): %s",
                  length(.unmapped), paste(sort(.unmapped), collapse = ", ")))
}
message(sprintf("Loaded %d rows, %d studies, %d features",
                nrow(df), length(unique(df$study_id)), length(unique(df$feature_id))))

# ---------------------------------------------------------------------------
# Per-model meta-analysis
# ---------------------------------------------------------------------------

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
    rma(yi = effect_size, sei = se, data = sub, method = method,
        slab = sub$study_id),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  data.frame(
    feature_id = sub$feature_id[1], k = res$k, study_ids = study_ids,
    yi = res$b[1], se = res$se,
    ci_lb = res$ci.lb, ci_ub = res$ci.ub,
    pval = res$pval, I2 = res$I2, tau2 = res$tau2,
    stringsAsFactors = FALSE)
}

all_model_results <- list()

for (model in SUBGROUPS) {
  sub_df <- df[!is.na(df$pain_model) & df$pain_model == model, ]
  if (nrow(sub_df) == 0) next
  message(sprintf("\n=== %s: %d studies, %d features ===",
                  model, length(unique(sub_df$study_id)),
                  length(unique(sub_df$feature_id))))

  feats <- split(sub_df, sub_df$feature_id)
  res_list <- lapply(feats, meta_one_feature)
  res_list <- Filter(Negate(is.null), res_list)
  pooled <- bind_rows(res_list)
  if (nrow(pooled) == 0) next

  pooled <- pooled %>%
    arrange(pval) %>%
    mutate(pain_model = model)
  pooled <- adjust_pooled_fdr(pooled, min_k = MIN_STUDIES)

  out_csv <- file.path(OUT_DIR, sprintf("%s_pooled.csv", model))
  write.csv(pooled, out_csv, row.names = FALSE)
  # Top hit is reported from the adjusted family, so k=1 pass-through rows
  # (padj = NA) can never be named as the leading signal for a stratum.
  top <- pooled[!is.na(pooled$padj), ]
  if (nrow(top) == 0) {
    message(sprintf("  → %d features, no feature reaches k>=%d",
                    nrow(pooled), MIN_STUDIES))
  } else {
    message(sprintf("  → %d features (%d pooled, k>=%d), top hit: %s (padj=%.2e)",
                    nrow(pooled), nrow(top), MIN_STUDIES, top$feature_id[1], top$padj[1]))
  }

  all_model_results[[model]] <- pooled
}

# ---------------------------------------------------------------------------
# Between-model heterogeneity test (Q_between via meta-regression)
# ---------------------------------------------------------------------------
# For each feature with data in >=2 models, run rma.mv with model as moderator.

message("\n=== Between-model heterogeneity test ===")

# Build per-study table restricted to the analysed subgroups (drops "multi",
# "other", "FLIT", "CFA_orofacial" and unmapped studies).
clean_df <- df[df$pain_model %in% SUBGROUPS, ]

# Restrict to features appearing in >=2 different models
feat_models <- clean_df %>%
  filter(!is.na(se), se > 0) %>%
  group_by(feature_id) %>%
  summarise(n_models = n_distinct(pain_model), k_total = n(), .groups = "drop") %>%
  filter(n_models >= 2, k_total >= 4)

message(sprintf("  Features in >=2 models with k>=4: %d", nrow(feat_models)))

q_results <- vector("list", nrow(feat_models))
for (i in seq_len(nrow(feat_models))) {
  feat <- feat_models$feature_id[i]
  sub <- usable_se(clean_df[clean_df$feature_id == feat, ])
  # droplevels: a feature seen in 2 models must not carry all 9 SUBGROUPS
  # levels, or metafor silently drops the aliased columns and the reported
  # degrees of freedom stop tracking the design.
  sub$model_f <- droplevels(factor(sub$pain_model, levels = SUBGROUPS))

  # Null model (intercept only, random by study)
  res0 <- tryCatch(
    rma(yi = effect_size, sei = se, data = sub, method = METHOD),
    error = function(e) NULL)
  # Model with pain_model as moderator
  res1 <- tryCatch(
    rma(yi = effect_size, sei = se, mods = ~ model_f, data = sub, method = METHOD),
    error = function(e) NULL)

  if (is.null(res0) || is.null(res1)) next

  # Q_between = Q_residual(null) - Q_residual(full) ≈ QM (moderator test)
  q_results[[i]] <- data.frame(
    feature_id = feat,
    k_total    = res0$k,
    n_models   = feat_models$n_models[i],
    QM         = res1$QM,            # moderator Wald test stat
    QM_df      = moderator_df(res1), # moderator df (NOT the coefficient count)
    QM_pval    = res1$QMp,           # p-value for between-model differences
    I2_total   = res0$I2,
    tau2_total = res0$tau2,
    stringsAsFactors = FALSE
  )
}

q_df <- bind_rows(Filter(Negate(is.null), q_results))
if (nrow(q_df) > 0) {
  q_df <- q_df %>%
    arrange(QM_pval) %>%
    mutate(padj_QM = p.adjust(QM_pval, method = "BH"))
  write.csv(q_df, file.path(OUT_DIR, "between_model_Q.csv"), row.names = FALSE)
  message(sprintf("  Saved between-model Q test for %d features", nrow(q_df)))
  message(sprintf("  Features with significant between-model diff (padj<0.05): %d",
                  sum(q_df$padj_QM < 0.05, na.rm = TRUE)))
  message("  Top 10 features with strongest between-model heterogeneity:")
  print(head(q_df[, c("feature_id","k_total","n_models","QM","QM_pval","padj_QM","I2_total")], 10))
}

# ---------------------------------------------------------------------------
# Known pain gene forest plots across models
# ---------------------------------------------------------------------------

message("\n=== Known gene per-model forest plots ===")

KNOWN_GENES <- c("Jun", "Cd44", "Atf3", "Tspo", "Aif1", "Cklf",
                 "Ccl2", "Bdnf", "Hpca", "RT1-DMb", "Rnaset2-ps1")

pdf_path <- file.path(OUT_DIR, "known_genes_forest.pdf")
pdf(pdf_path, width = 9, height = 6)

for (gene in KNOWN_GENES) {
  gene_df <- usable_se(clean_df[clean_df$feature_id == gene, ])
  if (nrow(gene_df) < 2) next

  # Per-model pooled estimates for the diamond
  model_ests <- do.call(rbind, lapply(SUBGROUPS, function(m) {
    sub <- gene_df[gene_df$pain_model == m, ]
    if (nrow(sub) < 2) return(NULL)
    res <- tryCatch(rma(yi = effect_size, sei = se, data = sub, method = METHOD),
                    error = function(e) NULL)
    if (is.null(res)) return(NULL)
    data.frame(label = sprintf("%s (k=%d)", m, res$k),
               yi = res$b[1], ci_lb = res$ci.lb, ci_ub = res$ci.ub,
               pval = res$pval, is_model = TRUE, stringsAsFactors = FALSE)
  }))

  # Individual studies
  ind_df <- gene_df %>%
    mutate(label    = study_id,
           yi       = effect_size,
           ci_lb    = effect_size - Z_975 * se,
           ci_ub    = effect_size + Z_975 * se,
           is_model = FALSE) %>%
    arrange(pain_model, effect_size) %>%
    select(label, yi, ci_lb, ci_ub, pval, is_model, pain_model)

  # model_ests is NULL when no single model reaches k>=2 for this gene (e.g. a
  # gene in two studies drawn from two different models, which still clears the
  # nrow(gene_df) >= 2 gate above). `NULL %>% mutate()` errors, and this line
  # sits outside the tryCatch below, so the script would abort with the pdf
  # device still open.
  all_rows <- if (is.null(model_ests)) {
    ind_df
  } else {
    bind_rows(model_ests %>% mutate(pain_model = NA), ind_df)
  }
  all_rows$y_pos <- seq_len(nrow(all_rows))

  tryCatch({
    par_default <- par(no.readonly = TRUE)
    plot_df <- all_rows[rev(seq_len(nrow(all_rows))), ]
    plot_df$y_rev <- seq_len(nrow(plot_df))

    xlim_val <- range(c(plot_df$ci_lb, plot_df$ci_ub), na.rm = TRUE)
    xlim_val[1] <- min(xlim_val[1], -0.5)
    xlim_val[2] <- max(xlim_val[2],  0.5)

    plot(NA, xlim = xlim_val, ylim = c(0.5, nrow(plot_df) + 0.5),
         xlab = "log2 Fold Change (case vs control)",
         yaxt = "n", ylab = "",
         main = sprintf("%s — per-model stratified forest plot", gene),
         cex.main = 0.9)
    abline(v = 0, lty = 2, col = "grey50")

    for (j in seq_len(nrow(plot_df))) {
      row <- plot_df[j, ]
      if (isTRUE(row$is_model)) {
        # Diamond for model-level pooled
        polygon(c(row$ci_lb, row$yi, row$ci_ub, row$yi),
                c(j, j + 0.25, j, j - 0.25),
                col = "#d6604d", border = "#d6604d")
        text(xlim_val[1], j, row$label, adj = 0, cex = 0.7, font = 2)
      } else {
        arrows(row$ci_lb, j, row$ci_ub, j, angle = 90, code = 3,
               length = 0.04, col = "steelblue", lwd = 0.8)
        points(row$yi, j, pch = 15, col = "steelblue", cex = 0.7)
        lbl <- sprintf("%s [%s]", row$label,
                       ifelse(is.na(row$pain_model), "", row$pain_model))
        text(xlim_val[1], j, lbl, adj = 0, cex = 0.6)
      }
    }
  }, error = function(e) message("  Forest error for ", gene, ": ", e$message))
}

dev.off()
message(sprintf("Known gene forest plots → %s", pdf_path))

# ---------------------------------------------------------------------------
# Cross-model effect heatmap for top multi-model features
# ---------------------------------------------------------------------------

message("\n=== Cross-model heatmap ===")

# Build a wide matrix: feature x model (yi values for features in >=2 models).
# When no feature qualified for the Q test, q_df is a 0-column tibble and
# filtering on padj_QM would error on the missing column.
top_features <- if (nrow(q_df) > 0 && "padj_QM" %in% colnames(q_df)) {
  q_df %>%
    filter(!is.na(padj_QM)) %>%
    arrange(QM_pval) %>%
    head(30) %>%
    pull(feature_id)
} else {
  character(0)
}

if (length(top_features) > 0) {
  wide <- do.call(rbind, lapply(all_model_results, function(res_df) {
    res_df[res_df$feature_id %in% top_features,
           c("feature_id", "pain_model", "yi")]
  }))
  if (!is.null(wide) && nrow(wide) > 0) {
    mat <- tidyr::pivot_wider(wide, names_from = "pain_model", values_from = "yi")
    mat_num <- as.matrix(mat[, -1])
    rownames(mat_num) <- mat$feature_id

    # Clamp for color scale
    mat_num[mat_num >  3] <-  3
    mat_num[mat_num < -3] <- -3

    # heatmap() clusters through dist(), which yields NaN for any pair of rows
    # sharing no observed model -- and hclust then aborts, taking the rest of
    # this script with it. With 14 strata, several of them human-only, that is
    # the normal case rather than the edge case.
    clusterable <- drop_sparse_for_clustering(mat_num)
    pdf(file.path(OUT_DIR, "model_heatmap.pdf"), width = 6, height = 8)
    heatmap_col <- colorRampPalette(c("#4575b4","white","#d73027"))(100)
    if (!is.null(clusterable)) {
      heatmap(clusterable, col = heatmap_col, scale = "none",
              main = "Top between-model features\n(log2FC, clamped ±3)",
              cexRow = 0.7, cexCol = 0.9, margins = c(5, 8))
      message(sprintf("Heatmap over %d features x %d models",
                      nrow(clusterable), ncol(clusterable)))
    } else {
      heatmap(mat_num, col = heatmap_col, scale = "none", Rowv = NA, Colv = NA,
              main = "Top between-model features\n(log2FC, clamped ±3; unclustered)",
              cexRow = 0.7, cexCol = 0.9, margins = c(5, 8))
      message("Matrix too sparse to cluster — heatmap drawn unclustered")
    }
    dev.off()
    message(sprintf("Heatmap → %s", file.path(OUT_DIR, "model_heatmap.pdf")))
  }
}

# ---------------------------------------------------------------------------
# Combined significant-feature table
# ---------------------------------------------------------------------------
# This file used to sit in results/ without any script producing it, which made
# it silently outlive the analyses it came from. It is now regenerated here on
# every run so its contents always match the per-stratum CSVs beside it.

sig_all <- bind_rows(lapply(all_model_results, function(res) {
  res[!is.na(res$padj) & res$padj < 0.05 & res$k >= 2, ]
}))
sig_all_path <- file.path(OUT_DIR, "significant_all_models.csv")
write.csv(sig_all, sig_all_path, row.names = FALSE)
message(sprintf("Significant features across all strata: %d → %s",
                nrow(sig_all), sig_all_path))

# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------

cat("\n========== Stratified Meta-Analysis Summary ==========\n")
for (model in names(all_model_results)) {
  res <- all_model_results[[model]]
  sig <- res[!is.na(res$padj) & res$padj < 0.05 & res$k >= MIN_STUDIES, ]
  cat(sprintf("\n[%s] %d features pooled, %d significant (k>=%d, padj<0.05)\n",
              model, nrow(res), MIN_STUDIES, nrow(sig)))
  if (nrow(sig) > 0) {
    print(head(sig[order(sig$padj), c("feature_id","k","yi","I2","padj")], 5))
  }
}

cat("\n[Between-model Q test] Top features with model-specific effects:\n")
if (exists("q_df") && nrow(q_df) > 0) {
  top_q <- q_df[!is.na(q_df$padj_QM) & q_df$padj_QM < 0.05, ]
  cat(sprintf("  %d features show significant between-model differences (padj<0.05)\n",
              nrow(top_q)))
  if (nrow(top_q) > 0) {
    print(head(top_q[, c("feature_id","k_total","n_models","QM_pval","padj_QM","I2_total")], 10))
  }
}

message("\nStep 06b complete. Output: ", OUT_DIR)
