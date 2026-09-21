#' Per-study differential analysis for proteomics (PRIDE Archive).
#'
#' Expects label-free quantification (LFQ) intensity matrices in tabular form.
#' Uses limma on log2-transformed LFQ intensities.
#'
#' Requires: limma, dplyr

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

# Directory of *this* file, for sourcing its siblings.
#
# These files live in R/ but are sourced by scripts in pipeline/, and the
# previous idiom derived the directory from commandArgs("--file="), which names
# the *invoking* script. Called as `Rscript pipeline/05_per_study_da.R` that
# resolved to pipeline/se_utils.R, which does not exist. It went unnoticed
# because each source() is guarded by an `exists()` check, so it only fires
# when nothing has already defined the function and the failure depends on
# source order.
#
# When a file is source()d, R records its path as `ofile` on the sourcing
# frame; walking the stack finds it. The commandArgs fallback still applies
# when this file is itself the script being run.
.sibling_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("--file=", "", fa[1]))))
  getwd()
}

if (!exists("assign_groups")) source(file.path(.sibling_dir(), "da_utils.R"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

#' Run limma DA on a proteomics LFQ matrix.
#'
#' @param accession   PRIDE accession (e.g. "PXD000001").
#' @param expr_matrix Matrix of log2-LFQ intensities (proteins × samples).
#'                    Rows are protein IDs (UniProt preferred), columns are samples.
#' @param sample_meta data.frame with a `group` column ("case" | "control").
#' @param cfg         Named list from conf/analysis/default.yaml.
#' @return data.frame of effect sizes or NULL.
run_proteomics_da <- function(accession, expr_matrix, sample_meta, cfg) {
  message(sprintf("\n=== %s [Proteomics] ===", accession))

  # Align samples
  common_samples <- intersect(colnames(expr_matrix), rownames(sample_meta))
  if (length(common_samples) == 0) {
    message("  No common samples between matrix and metadata.")
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

  usable  <- pheno[!is.na(pheno$group), ]
  expr_u  <- expr[, rownames(usable), drop = FALSE]
  group   <- factor(usable$group, levels = c("control", "case"))
  design  <- model.matrix(~group)

  # Remove proteins with > 50% missing across samples
  na_frac <- rowMeans(is.na(expr_u))
  expr_u  <- expr_u[na_frac <= 0.5, , drop = FALSE]

  # Impute remaining NAs with half-minimum (common in proteomics)
  expr_u <- apply(expr_u, 2, function(col) {
    col[is.na(col)] <- min(col, na.rm = TRUE) / 2
    col
  })

  fit  <- lmFit(expr_u, design)
  fit2 <- eBayes(fit)
  tt   <- topTable(fit2, coef = "groupcase", number = Inf, sort.by = "none")

  se_approx <- abs(tt$logFC / tt$t)

  build_effect_df(rownames(tt), tt$logFC, se_approx,
                  tt$P.Value, tt$adj.P.Val,
                  n_case, n_ctrl, accession,
                  "proteomics", "UniProt", "log2FC")
}


# ---------------------------------------------------------------------------
# Loader helper (for PRIDE tabular exports)
# ---------------------------------------------------------------------------

#' Load a PRIDE expression matrix from a tabular file (TSV/CSV).
#'
#' Expects the first column to contain protein IDs; remaining columns are samples.
#' Intensities are log2-transformed if they appear to be raw (max > 100).
#'
#' @param file_path Path to the expression file.
#' @return Numeric matrix (proteins × samples) or NULL.
load_pride_matrix <- function(file_path) {
  if (!file.exists(file_path)) {
    message("  Expression file not found: ", file_path)
    return(NULL)
  }
  sep <- if (grepl("\\.tsv$", file_path)) "\t" else ","
  df  <- read.table(file_path, sep = sep, header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)

  id_col  <- df[[1]]
  mat     <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- id_col

  # Auto log2 transform
  if (max(mat, na.rm = TRUE) > 100) {
    mat <- log2(mat + 1)
    message("  Applied log2 transformation (intensities were raw).")
  }
  mat
}
