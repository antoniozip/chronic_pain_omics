#' Per-study differential analysis for metabolomics and lipidomics.
#'
#' Uses limma on log2-transformed concentration data.
#' Standardized mean difference (SMD / Hedges' g) is computed for
#' inter-study comparability.
#'
#' Requires: limma, dplyr

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

if (!exists("assign_groups")) source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))), "da_utils.R"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

#' Run differential analysis for one metabolomics or lipidomics study.
#'
#' @param accession   MetaboLights accession (e.g. "MTBLS123").
#' @param expr_matrix Matrix of metabolite/lipid concentrations (features × samples).
#'                    Values are assumed to be on a normalized concentration scale.
#' @param sample_meta data.frame with a `group` column ("case" | "control").
#' @param id_space    "ChEBI", "HMDB", or "LIPIDMAPS".
#' @param modality    "metabolomics" or "lipidomics".
#' @param cfg         Named list from conf/analysis/default.yaml.
#' @return data.frame of effect sizes or NULL.
run_metabolomics_da <- function(accession, expr_matrix, sample_meta,
                                id_space, modality = "metabolomics", cfg) {
  message(sprintf("\n=== %s [%s] ===", accession, modality))

  common_samples <- intersect(colnames(expr_matrix), rownames(sample_meta))
  if (length(common_samples) == 0) {
    message("  No common samples.")
    return(NULL)
  }

  expr  <- expr_matrix[, common_samples, drop = FALSE]
  pheno <- sample_meta[common_samples, , drop = FALSE]
  pheno <- assign_groups(pheno, cfg$da$case_labels, cfg$da$control_labels)

  n_case <- sum(pheno$group == "case", na.rm = TRUE)
  n_ctrl <- sum(pheno$group == "control", na.rm = TRUE)

  if (n_case < cfg$da$min_samples_per_group || n_ctrl < cfg$da$min_samples_per_group) {
    message(sprintf("  Insufficient group sizes (case=%d, ctrl=%d)", n_case, n_ctrl))
    return(NULL)
  }

  usable <- pheno[!is.na(pheno$group), , drop = FALSE]
  expr_u <- expr[, rownames(usable), drop = FALSE]
  group  <- factor(usable$group, levels = c("control", "case"))

  # Log2 transform if data appears raw
  if (max(expr_u, na.rm = TRUE) > 100) {
    expr_u <- log2(expr_u + 1)
  }

  # Remove features with > 30% missing
  na_frac <- rowMeans(is.na(expr_u))
  expr_u  <- expr_u[na_frac <= 0.3, , drop = FALSE]

  design <- model.matrix(~group)
  fit    <- lmFit(expr_u, design)
  fit2   <- eBayes(fit, trend = TRUE)
  tt     <- topTable(fit2, coef = "groupcase", number = Inf, sort.by = "none")

  # Compute Hedges' g (SMD) for cross-study comparison
  smd <- compute_hedges_g(expr_u, group)
  smd_matched <- smd[match(rownames(tt), names(smd))]

  se_approx <- abs(tt$logFC / tt$t)

  df <- build_effect_df(rownames(tt), smd_matched, se_approx,
                        tt$P.Value, tt$adj.P.Val,
                        n_case, n_ctrl, accession,
                        modality, id_space, "SMD")

  # Also store raw log2FC for reference
  df$log2fc_raw <- tt$logFC
  df
}


# ---------------------------------------------------------------------------
# Hedges' g
# ---------------------------------------------------------------------------

#' Compute Hedges' g for each feature (row) given case/control group vector.
#'
#' @param mat   Numeric matrix (features × samples).
#' @param group Factor with levels c("control", "case").
#' @return Named numeric vector of Hedges' g values.
compute_hedges_g <- function(mat, group) {
  n1 <- sum(group == "case")
  n0 <- sum(group == "control")
  # Correction factor J (Hedges 1981)
  df_total <- n1 + n0 - 2
  j <- 1 - (3 / (4 * df_total - 1))

  apply(mat, 1, function(row) {
    m1 <- mean(row[group == "case"],    na.rm = TRUE)
    m0 <- mean(row[group == "control"], na.rm = TRUE)
    s1 <- var(row[group == "case"],    na.rm = TRUE)
    s0 <- var(row[group == "control"], na.rm = TRUE)
    sp <- sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / df_total)
    if (is.na(sp) || sp == 0) return(NA_real_)
    j * (m1 - m0) / sp
  })
}


# ---------------------------------------------------------------------------
# MetaboLights loader helper
# ---------------------------------------------------------------------------

#' Load a MetaboLights data matrix (m_*.tsv) and sample metadata (s_*.tsv).
#'
#' @param study_dir Path to the downloaded MetaboLights study directory.
#' @return List with `matrix` (features × samples) and `sample_meta` (data.frame).
load_metabolights_study <- function(study_dir) {
  data_files   <- list.files(study_dir, pattern = "^m_.*\\.tsv$", full.names = TRUE)
  sample_files <- list.files(study_dir, pattern = "^s_.*\\.tsv$", full.names = TRUE)

  if (length(data_files) == 0 || length(sample_files) == 0) {
    message("  Missing m_*.tsv or s_*.tsv in ", study_dir)
    return(NULL)
  }

  mat_df <- read.table(data_files[1], sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- mat_df[[1]]
  mat    <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- id_col

  sample_df <- read.table(sample_files[1], sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sample_df) <- sample_df[[1]]

  list(matrix = mat, sample_meta = sample_df)
}
