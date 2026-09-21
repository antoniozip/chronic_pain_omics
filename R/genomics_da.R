#' Per-study effect-size extraction for genomics (GWAS Catalog).
#'
#' GWAS Catalog studies already provide summary statistics via the REST API.
#' This module fetches associations for each included study and standardizes
#' them into the shared effect-size format for meta-analysis.
#'
#' Requires: httr, jsonlite, dplyr

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
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

GWAS_API <- "https://www.ebi.ac.uk/gwas/rest/api"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

#' Extract GWAS summary statistics for one GWAS Catalog study.
#'
#' @param accession   GWAS Catalog accession (e.g. "GCST000001").
#' @param effect_unit Inferred unit from harmonized record ("beta" | "OR").
#' @param cfg         Named list from conf/analysis/default.yaml.
#' @return data.frame of effect sizes or NULL.
run_genomics_da <- function(accession, effect_unit, cfg) {
  message(sprintf("\n=== %s [GWAS] ===", accession))

  associations <- fetch_gwas_associations(accession)
  if (is.null(associations) || nrow(associations) == 0) {
    message("  No associations retrieved.")
    return(NULL)
  }

  study_meta <- fetch_gwas_study_meta(accession)
  n_case    <- parse_sample_count(study_meta$initialSampleSize)
  n_control <- parse_sample_count(study_meta$replicationSampleSize)

  df <- standardize_gwas(associations, effect_unit)
  if (nrow(df) == 0) return(NULL)

  build_effect_df(df$feature_id, df$effect_size, df$se,
                  df$pval, df$pval,  # GWAS Catalog has no padj; use pval twice
                  n_case, n_control, accession,
                  "genomics", "rsID", effect_unit)
}


# ---------------------------------------------------------------------------
# GWAS Catalog REST helpers
# ---------------------------------------------------------------------------

fetch_gwas_associations <- function(accession, page_size = 500L) {
  url <- sprintf("%s/studies/%s/associations", GWAS_API, accession)
  resp <- tryCatch(
    GET(url, query = list(size = page_size, projection = "associationByStudy")),
    error = function(e) { message("  GWAS API error: ", e$message); NULL }
  )
  if (is.null(resp) || http_error(resp)) return(NULL)

  content_list <- content(resp, as = "parsed", type = "application/json")
  assoc_list   <- content_list[["_embedded"]][["associations"]]
  if (is.null(assoc_list)) return(NULL)

  bind_rows(lapply(assoc_list, parse_association))
}


fetch_gwas_study_meta <- function(accession) {
  url  <- sprintf("%s/studies/%s", GWAS_API, accession)
  resp <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(resp) || http_error(resp)) return(list())
  content(resp, as = "parsed", type = "application/json")
}


parse_association <- function(a) {
  snps <- a[["loci"]][[1]][["strongestRiskAlleles"]]
  rsid <- if (!is.null(snps)) snps[[1]][["riskAlleleName"]] else NA_character_
  rsid <- sub("-.*", "", rsid %||% NA_character_)  # strip allele suffix

  data.frame(
    rsid         = rsid,
    beta_or      = a[["betaNum"]]  %||% NA_real_,
    beta_unit    = a[["betaUnit"]] %||% NA_character_,
    or_per_copy  = a[["orPerCopyNum"]] %||% NA_real_,
    pval         = a[["pvalueMantissa"]] * 10^a[["pvalueExponent"]],
    stringsAsFactors = FALSE
  )
}


standardize_gwas <- function(df, effect_unit) {
  df <- df[!is.na(df$rsid), ]
  if (nrow(df) == 0) return(df)

  if (effect_unit == "OR") {
    df$effect_size <- log(df$or_per_copy)          # log(OR) for meta-analysis
    df$se <- NA_real_                               # GWAS Catalog lacks SE; imputed from CI if available
  } else {
    df$effect_size <- df$beta_or
    df$se <- NA_real_
  }

  # Flag rows with missing SE — meta-analysis will use p-value-derived SE if needed
  df$feature_id <- df$rsid
  df[c("feature_id", "effect_size", "se", "pval")]
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

parse_sample_count <- function(descriptor) {
  if (is.null(descriptor) || !nzchar(descriptor)) return(0L)
  m <- regmatches(descriptor, regexpr("[0-9,]+", descriptor))
  if (length(m) == 0) return(0L)
  as.integer(gsub(",", "", m[1]))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
