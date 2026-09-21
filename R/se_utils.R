#' Pure standard-error helpers shared by the differential-analysis backends.
#'
#' These formulas are the integration currency of the meta-analysis: a wrong SE
#' silently corrupts every pooled estimate downstream. They live in their own
#' file, free of Bioconductor dependencies, so that they remain unit-testable.
#'
#' Requires: nothing.

#' SE from a moderated t-statistic (limma).
#'
#' SE = |logFC / t|. NOT logFC / qt(0.975, df).
#'
#' @param logFC Numeric vector of log2 fold changes.
#' @param t     Numeric vector of t-statistics.
#' @return Numeric vector of standard errors.
se_from_t <- function(logFC, t) {
  abs(logFC / t)
}

#' SE back-calculated from a two-sided p-value (Cuffdiff).
#'
#' Cuffdiff's `test_stat` is NOT logFC/SE, so it cannot be used directly.
#' Invert the Normal approximation instead: SE = |effect| / qnorm(1 - p/2).
#' p is clamped away from 0 and 1 to avoid Inf and 0.
#'
#' @param effect Numeric vector of effect sizes (log2 fold changes).
#' @param pval   Numeric vector of two-sided p-values.
#' @return Numeric vector of standard errors; NA where undefined.
se_from_pvalue <- function(effect, pval) {
  p_clamped <- pmax(pmin(pval, 1 - 1e-10), 1e-10)
  z_abs     <- qnorm(1 - p_clamped / 2)
  ifelse(
    z_abs > 0.01 & !is.na(effect) & abs(effect) > 0,
    abs(effect) / z_abs,
    NA_real_
  )
}

#' Detect linear-scale (non-log2) expression data.
#'
#' GEO series matrices sometimes store raw or MAS5 intensities rather than log2
#' values. Threshold is Q75 > 50, deliberately NOT 100: GSE10238 has Q75 = 92.2
#' and must be transformed.
#'
#' @param expr      Numeric matrix or vector of expression values.
#' @param threshold Q75 above which data is assumed linear-scale.
#' @return TRUE if the data should be log2-transformed.
needs_log2 <- function(expr, threshold = 50) {
  q75 <- quantile(expr, 0.75, na.rm = TRUE)
  unname(q75 > threshold)
}

#' Apply log2(pmax(x, 1)) iff the data appears to be on a linear scale.
#'
#' @param expr      Numeric matrix or vector of expression values.
#' @param threshold Passed to needs_log2().
#' @return expr, log2-transformed if required.
maybe_log2 <- function(expr, threshold = 50) {
  if (needs_log2(expr, threshold)) {
    return(log2(pmax(expr, 1)))
  }
  expr
}
