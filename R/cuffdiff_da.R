#' Cuffdiff differential-expression output parsing.
#'
#' Kept separate from transcriptomics_da.R because it needs no Bioconductor
#' packages, and is therefore unit-testable in environments without them.
#'
#' Requires: da_utils.R (build_effect_df), se_utils.R (se_from_pvalue)

#' Parse a Cuffdiff gene_exp.diff file into a standardized effect data.frame.
#'
#' @param file_path Path to the Cuffdiff differential-expression table.
#' @param accession Study accession, written into study_id.
#' @param cfg       Config list (unused; kept for backend signature symmetry).
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

  d <- d[d$status %in% c("OK", "NOTEST"), ]
  d <- d[d$gene != "-" & nchar(d$gene) > 0, ]
  if (nrow(d) == 0) return(NULL)

  fc   <- as.numeric(d[["log2(fold_change)"]])
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

  # Cuffdiff emits one row per isoform/locus; keep the most significant per gene.
  tmp <- data.frame(
    gene = d$gene, fc = fc, se = se_approx,
    pval = pval, qval = as.numeric(d$q_value),
    stringsAsFactors = FALSE
  )
  tmp <- tmp[order(tmp$pval), ]
  tmp <- tmp[!duplicated(tmp$gene), ]

  message(sprintf("  Cuffdiff DE output: %d genes (after dedup)", nrow(tmp)))
  build_effect_df(
    tmp$gene, tmp$fc, tmp$se,
    tmp$pval, tmp$qval,
    NA_integer_, NA_integer_,
    accession, "transcriptomics", "symbol", "log2FC"
  )
}
