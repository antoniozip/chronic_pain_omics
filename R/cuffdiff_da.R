#' Cuffdiff differential-expression output parsing.
#'
#' Kept separate from transcriptomics_da.R because it needs no Bioconductor
#' packages, and is therefore unit-testable in environments without them.
#'
#' Requires: da_utils.R (build_effect_df), se_utils.R (se_from_pvalue)

#' Which way round the Cuffdiff contrast runs.
#'
#' Cuffdiff's `log2(fold_change)` is `log2(value_2 / value_1)`, i.e. sample_2
#' over sample_1 — the names are in the same file and are the only evidence of
#' which arm was the denominator. Every other backend fixes case-as-positive
#' (limma `coef = "groupcase"`, metabolomics `m1 - m0`, single cell
#' `contrast = c("group", "case", "control")`), so a control-vs-case invocation
#' silently entered the pool with case in the denominator and a sign-inverted
#' effect that cancels the other studies on a shared feature rather than
#' reinforcing them. `se_from_pvalue()` is symmetric and survives either sign,
#' so nothing else notices.
#'
#' @param sample_1,sample_2 Cuffdiff's two arm-name columns.
#' @param cfg Config list; `cfg$da$case_labels`/`control_labels` extend the
#'   literal "case"/"control" vocabulary.
#' @return Numeric vector: +1 where sample_2 is the case arm, -1 where it is
#'   the control arm. Stops when a row names neither, or both.
cuffdiff_contrast_direction <- function(sample_1, sample_2, cfg) {
  labels <- function(which, literal) {
    extra <- if (!is.null(cfg$da)) cfg$da[[which]] else NULL
    tolower(c(literal, extra))
  }
  hit <- function(x, words) {
    x <- tolower(as.character(x))
    out <- rep(FALSE, length(x))
    for (w in words) out <- out | grepl(w, x, fixed = TRUE)
    out
  }

  case_words <- labels("case_labels", "case")
  ctrl_words <- labels("control_labels", "control")

  s1_case <- hit(sample_1, case_words); s1_ctrl <- hit(sample_1, ctrl_words)
  s2_case <- hit(sample_2, case_words); s2_ctrl <- hit(sample_2, ctrl_words)

  keep <- s2_case & s1_ctrl          # sample_2 is case: sign is already right
  flip <- s2_ctrl & s1_case          # sample_2 is control: sign is inverted
  bad  <- (!keep & !flip) |
    (s1_case & s1_ctrl) | (s2_case & s2_ctrl)   # neither arm, or both

  if (any(bad)) {
    pairs <- unique(sprintf("%s vs %s", sample_1[bad], sample_2[bad]))
    stop("Cuffdiff contrast direction is not established for: ",
         paste(pairs, collapse = "; "),
         ". The project's convention is case minus control; refusing to guess ",
         "the sign. Rename the arms to name a case and a control.")
  }
  ifelse(flip, -1, 1)
}


#' @param file_path Path to the Cuffdiff differential-expression table.
#' @param accession Study accession, written into study_id.
#' @param cfg       Config list; its case/control label vocabulary resolves the
#'                  contrast direction (see cuffdiff_contrast_direction()).
#' @return data.frame from build_effect_df(), or NULL if the file is unusable.
run_cuffdiff_output <- function(file_path, accession, cfg) {
  d <- tryCatch(
    read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      message("  Cuffdiff read error: ", e$message)
      NULL
    }
  )
  if (is.null(d)) return(NULL)

  needed <- c("gene", "log2(fold_change)", "p_value", "q_value", "test_stat", "status")
  if (!all(needed %in% colnames(d))) return(NULL)

  # The arm names are the only evidence of which side was the denominator, so
  # their absence is not a "file we cannot use" case that may pass silently:
  # dropping the study quietly is how it would have entered with an arbitrary
  # sign before.
  if (!all(c("sample_1", "sample_2") %in% colnames(d))) {
    stop("Cuffdiff table has no sample_1/sample_2 columns, so the contrast ",
         "direction (case minus control) cannot be established: ", file_path)
  }

  d <- d[d$status %in% c("OK", "NOTEST"), ]
  d <- d[d$gene != "-" & nchar(d$gene) > 0, ]
  if (nrow(d) == 0) return(NULL)

  # Establish case-as-positive before anything is signed. A file whose arms
  # name neither case nor control stops here rather than entering the pool
  # with an arbitrary sign.
  direction <- cuffdiff_contrast_direction(d$sample_1, d$sample_2, cfg)

  fc   <- direction * as.numeric(d[["log2(fold_change)"]])
  pval <- as.numeric(d[["p_value"]])

  # Cuffdiff's test_stat is NOT logFC/SE. Invert the p-value instead.
  se_approx <- se_from_pvalue(fc, pval)

  ok <- !is.na(fc) & !is.na(se_approx) & is.finite(fc) &
        is.finite(se_approx) & se_approx > 0
  d         <- d[ok, ]
  fc        <- fc[ok]
  pval      <- pval[ok]
  se_approx <- se_approx[ok]
  if (nrow(d) == 0) return(NULL)

  # Cuffdiff emits one row per isoform/locus. Collapse to one row per gene by
  # largest |log2FC| — a pre-specified rule that does not look at the p-value.
  # The previous rule was smallest p-value, i.e. selection on the outcome: 05b
  # removed the identical rule from aggregate_to_gene() because the minimum of
  # n p-values is not uniform under the null and the retained effect is biased
  # away from zero. Ties break on the gene's row order, so the choice is
  # reproducible.
  tmp <- data.frame(
    gene = d$gene, fc = fc, se = se_approx,
    pval = pval, qval = as.numeric(d$q_value),
    stringsAsFactors = FALSE
  )
  tmp <- tmp[order(-abs(tmp$fc)), ]
  tmp <- tmp[!duplicated(tmp$gene), ]

  message(sprintf("  Cuffdiff DE output: %d genes (after dedup)", nrow(tmp)))
  build_effect_df(
    tmp$gene, tmp$fc, tmp$se,
    tmp$pval, tmp$qval,
    NA_integer_, NA_integer_,
    accession, "transcriptomics", "symbol", "log2FC"
  )
}
