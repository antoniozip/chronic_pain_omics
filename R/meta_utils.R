#' Shared utilities for the meta-analysis steps (06, 06b, 06c, 06d).
#'
#' Multiple-testing correction is the piece most easily got wrong here, so it
#' lives in one tested place rather than being re-implemented per script.


#' qnorm(0.975): converts a 95% CI half-width to a standard error, and back.
#'
#' Mirrors `_Z_975` in src/cp_multiomics/genomics/harmonize.py so the R and
#' Python sides of the pipeline place interval bounds identically.
Z_975 <- 1.959964


# ---------------------------------------------------------------------------
# Study selection
# ---------------------------------------------------------------------------

#' Accessions superseded by a finer-grained re-analysis of the same samples.
#'
#' When a study is split (e.g. by tissue) the original combined table must not
#' also enter the meta-analysis, or those samples are pooled twice and their
#' correlated estimates are treated as independent. The registry records the
#' replacement and the reason for each exclusion, so no study is dropped
#' without a documented decision (CLAUDE.md: never silently drop a study).
#'
#' @param path Path to conf/analysis/superseded_studies.csv.
#' @return Character vector of excluded study_ids; empty if the file is absent.
superseded_studies <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (nrow(df) == 0 || !"study_id" %in% colnames(df)) {
    return(character(0))
  }
  as.character(df$study_id)
}


# ---------------------------------------------------------------------------
# Row filtering
# ---------------------------------------------------------------------------

#' Keep only rows carrying a usable standard error.
#'
#' A meta-analysis needs a positive, non-missing SE for every contributing
#' estimate; rows failing that cannot be weighted and are dropped everywhere
#' this filter is applied.
#'
#' @param df     data.frame of per-study or per-feature estimates.
#' @param se_col Name of the standard-error column.
#' @return `df` restricted to rows where `se_col` is non-missing and > 0.
usable_se <- function(df, se_col = "se") {
  if (!se_col %in% colnames(df)) {
    stop("usable_se: missing column: ", se_col)
  }
  se <- df[[se_col]]
  df[!is.na(se) & se > 0, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Feature selection for reporting
# ---------------------------------------------------------------------------

#' Pick the strongest features that were actually meta-analysed.
#'
#' Ranking over every row selects single-study pass-throughs, whose raw
#' p-values carry no heterogeneity penalty and are routinely smaller than any
#' pooled p-value. A forest plot cannot draw those, so they are silently
#' skipped downstream and the figure renders fewer panels than requested.
#' Selecting within the adjusted family avoids both problems.
#'
#' @param df       data.frame of pooled per-feature results.
#' @param n        Maximum number of features to return.
#' @param pval_col Name of the raw pooled p-value column.
#' @param padj_col Adjusted p-value column; non-NA marks the pooled family.
#' @return Character vector of feature ids, ascending by p-value, length <= n.
top_pooled_features <- function(df, n, pval_col = "pval_pooled",
                                padj_col = "padj_pooled") {
  missing_cols <- setdiff(c("feature_id", pval_col, padj_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop("top_pooled_features: missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  pooled <- df[!is.na(df[[padj_col]]), , drop = FALSE]
  if (nrow(pooled) == 0) return(character(0))
  ranked <- pooled$feature_id[order(pooled[[pval_col]])]
  utils::head(as.character(ranked), n)
}


# ---------------------------------------------------------------------------
# Meta-regression reporting
# ---------------------------------------------------------------------------

#' Moderator degrees of freedom for a metafor meta-regression.
#'
#' `res$p` counts the model's coefficients, intercept included, and is the
#' wrong quantity to publish as a moderator df: with g groups it reports g
#' rather than g-1. metafor exposes the correct value as `QMdf[1]`, which also
#' stays correct when rank-deficient columns are dropped from the design.
#'
#' @param res An rma object fitted with `mods`.
#' @return Moderator degrees of freedom for the QM test.
moderator_df <- function(res) {
  unname(res$QMdf[1])
}


# ---------------------------------------------------------------------------
# Multiple-testing correction
# ---------------------------------------------------------------------------

#' Adjust p-values within the pooled (k >= min_k) family only.
#'
#' Per-feature meta-analysis frames carry two kinds of row. Features seen in a
#' single study are passed through with that study's raw p-value; features seen
#' in >= 2 studies carry a genuine meta-analytic p-value. Only the latter are
#' the inferential family, and downstream reporting filters on `k >= 2`
#' anyway — so adjusting over both kinds both inflates m and lets single-study
#' p-values occupy the low ranks, shifting the adjusted p-value of every pooled
#' feature in an uncontrolled direction.
#'
#' Non-qualifying rows are kept with `out_col` = NA rather than dropped, so no
#' feature disappears from the output (CLAUDE.md: never silently drop a study).
#'
#' @param df        data.frame of per-feature pooled results.
#' @param pval_col  Name of the raw p-value column.
#' @param k_col     Name of the study-count column.
#' @param min_k     Minimum studies for a row to enter the family.
#' @param out_col   Name of the adjusted p-value column to add.
#' @param method    Any method accepted by stats::p.adjust.
#' @return `df` with `out_col` added; row order and row count are preserved.
adjust_pooled_fdr <- function(df, pval_col = "pval", k_col = "k", min_k = 2L,
                              out_col = "padj", method = "BH") {
  missing_cols <- setdiff(c(pval_col, k_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop("adjust_pooled_fdr: missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  pvals <- as.numeric(df[[pval_col]])
  ks    <- as.numeric(df[[k_col]])
  in_family <- !is.na(pvals) & !is.na(ks) & ks >= min_k

  df[[out_col]] <- rep(NA_real_, nrow(df))
  if (any(in_family)) {
    df[[out_col]][in_family] <- p.adjust(pvals[in_family], method = method)
  }
  df
}


#' Prepare a feature-by-model matrix for hierarchical clustering, or decline.
#'
#' `heatmap()` clusters rows and columns through `dist()`, which returns NaN
#' for any pair sharing no observed column, and `hclust()` then aborts with
#' "NA/NaN/Inf in foreign function call (arg 10)". That is not a hypothetical:
#' it ended the 2026-08-29 `06b` run after every stratum had been pooled and
#' written, because the search amendment took the stratified design from 9
#' models to 14 and a top between-model feature is now absent from most of the
#' new, human-only columns.
#'
#' Rows and columns carrying fewer than `min_obs` observations are dropped;
#' if what survives still contains a pair with no shared column, NULL is
#' returned so the caller can draw the heatmap unclustered rather than hand
#' `hclust` a matrix that will abort.
#'
#' @param mat     Numeric matrix, features in rows and models in columns.
#' @param min_obs Minimum non-NA entries a row or column must carry.
#' @return The trimmed matrix, or NULL if it cannot be clustered.
drop_sparse_for_clustering <- function(mat, min_obs = 2) {
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  keep_col <- colSums(!is.na(mat)) >= min_obs
  mat <- mat[, keep_col, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)
  keep_row <- rowSums(!is.na(mat)) >= min_obs
  mat <- mat[keep_row, , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)

  # Every pair must share at least one observed column, or dist() yields NaN.
  observed <- !is.na(mat)
  shared_rows <- observed %*% t(observed)
  if (any(shared_rows == 0)) return(NULL)
  shared_cols <- t(observed) %*% observed
  if (any(shared_cols == 0)) return(NULL)
  mat
}


# ---------------------------------------------------------------------------
# Derived study units
# ---------------------------------------------------------------------------

#' Look a study up in an accession-keyed map, falling back to its parent.
#'
#' A derived unit (`GSE241361_DRG`, `GSE180627_S1`) is in no accession-keyed
#' map by construction: only the parent accession is ever curated. Every such
#' map is therefore a silent hole once a study is split -- the lookup returns
#' NA and the unit leaves the analysis without a message. That is what
#' splitting GSE180627 and GSE197233 on 2026-08-30 did to `06b`'s STUDY_MODEL:
#' the two parents left the pool, their six regional units matched nothing, and
#' CCI went from 8 studies to 7 and SNI from 18 to 17 while gaining six units.
#'
#' Suffixes are stripped one `_<token>` at a time, so a unit derived twice
#' still resolves, and an exact hit always wins over an inherited one -- a
#' split that changes the stratum can still be recorded per unit.
#'
#' Inheriting is only right where the split does not change the mapped
#' property. Region does not change the pain model; it does change the tissue,
#' which is why `data/interim/sni_tissue_annotations.csv` names derived units
#' explicitly instead of relying on this.
#'
#' @param ids Character vector of study ids.
#' @param map Named character vector keyed by study id.
#' @return Character vector the length of `ids`; NA where nothing resolved.
resolve_by_parent <- function(ids, map) {
  ids <- as.character(ids)
  out <- unname(map[ids])
  if (length(ids) == 0) return(character(0))
  candidate <- ids
  repeat {
    miss <- which(is.na(out))
    if (length(miss) == 0) break
    parent <- sub("_[^_]+$", "", candidate[miss])
    if (all(parent == candidate[miss])) break
    candidate[miss] <- parent
    out[miss] <- unname(map[parent])
  }
  out
}
