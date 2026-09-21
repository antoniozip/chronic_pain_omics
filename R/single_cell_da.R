#' Per-study differential analysis for the pseudobulked single-cell arm.
#'
#' The whole design of the arm is that this file has nothing new in it. Step 19
#' sums each sample's cells into one counts column, which makes the sample the
#' unit of replication and the table a bulk table; from there the analysis is
#' the same DESeq2 Wald test the RNA-seq branch of `transcriptomics_da.R` runs,
#' on the same unshrunken standard error, emitting the same effect table
#' through the same `build_effect_df()`. Anything that diverged here would make
#' the arm's effect sizes incomparable with the pooled bulk estimates, which is
#' the one comparison the arm exists to support.
#'
#' Two things do differ, and both are inputs rather than method:
#'
#'   * Counts come from `data/interim/single_cell/pseudobulk/`, not from GEO.
#'     Nothing is downloaded here; step 19 already resolved four deposit
#'     formats into one.
#'   * Groups come from `data/interim/single_cell/groups/`, written by step 18
#'     from curated per-study rules. They are deliberately not the geo_cache
#'     `{accession}_groups.csv` the bulk arm reads: every one of these
#'     accessions is also in the bulk arm's screening sheet, excluded, and a
#'     groups file dropped into the shared cache is exactly the kind of thing
#'     that later gets picked up by the wrong pipeline.
#'
#' Requires: DESeq2.

suppressPackageStartupMessages({
  library(DESeq2)
})

.sc_sibling_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("--file=", "", fa[1]))))
  getwd()
}

if (!exists("build_effect_df")) source(file.path(.sc_sibling_dir(), "da_utils.R"))


#' Run differential analysis for one pseudobulked single-cell study.
#'
#' @param accession   GEO series accession.
#' @param cfg         Named list from conf/analysis/default.yaml.
#' @param interim_dir data/interim, holding single_cell/pseudobulk and groups.
#' @return data.frame of effect sizes, or NULL when the study cannot be run.
run_single_cell_da <- function(accession, cfg,
                               interim_dir = file.path("data", "interim")) {
  message(sprintf("\n=== %s [pseudobulk] ===", accession))

  base <- file.path(interim_dir, "single_cell")
  counts_path <- file.path(base, "pseudobulk", paste0(accession, "_counts.tsv.gz"))
  groups_path <- file.path(base, "groups", paste0(accession, "_groups.csv"))
  if (!file.exists(counts_path)) {
    message("  No pseudobulk table — run step 19 first"); return(NULL)
  }
  if (!file.exists(groups_path)) {
    message("  No group assignment — run step 18 first"); return(NULL)
  }

  counts <- read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE)
  groups <- read.csv(groups_path, stringsAsFactors = FALSE)

  # Step 19 sums the libraries of a repeated animal into one column named for
  # the subject, so the arm assignment has to be read at that level too.
  # Taking it per GSM would enter the animal twice under two different names,
  # which is the pseudoreplication the subject id exists to prevent.
  by_subject <- unique(groups[, c("subject", "group")])
  arm <- setNames(by_subject$group, by_subject$subject)
  present <- intersect(colnames(counts), names(arm))
  missing <- setdiff(names(arm), colnames(counts))
  if (length(missing) > 0) {
    # Not an error: a sample can be described in the series metadata and absent
    # from the deposit. Step 19 records which, and this is where the arm counts
    # actually shrink, so say so rather than let n drift silently.
    message(sprintf("  %d assigned sample(s) absent from the pseudobulk table: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }
  counts <- counts[, present, drop = FALSE]
  group_vec <- unname(arm[present])

  n_case <- sum(group_vec == "case")
  n_ctrl <- sum(group_vec == "control")
  message(sprintf("  %d case / %d control", n_case, n_ctrl))
  if (n_case < cfg$da$min_samples_per_group || n_ctrl < cfg$da$min_samples_per_group) {
    message("  Below da.min_samples_per_group — skipping"); return(NULL)
  }

  mat <- as.matrix(counts)
  mode(mat) <- "integer"
  if (anyNA(mat)) {
    message("  Non-integer or missing counts in the pseudobulk table"); return(NULL)
  }
  mat <- mat[rowSums(mat) >= cfg$da$min_count, , drop = FALSE]
  if (nrow(mat) == 0) { message("  No features pass min_count"); return(NULL) }

  col_data <- data.frame(group = factor(group_vec, levels = c("control", "case")),
                         row.names = colnames(mat))
  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = mat, colData = col_data, design = ~group),
    error = function(e) { message("  DESeq2 setup: ", e$message); NULL }
  )
  if (is.null(dds)) return(NULL)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "case", "control"),
                 alpha = cfg$da$padj_threshold)

  # Unshrunken MLE, as in the bulk arm: apeglm's posterior SD is two orders of
  # magnitude below the true SE and would dominate inverse-variance weighting.
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$lfcSE) &
                   res_df$lfcSE > 0, , drop = FALSE]
  message(sprintf("  DESeq2 pseudobulk: %d features", nrow(res_df)))
  if (nrow(res_df) == 0) return(NULL)

  build_effect_df(rownames(res_df), res_df$log2FoldChange, res_df$lfcSE,
                  res_df$pvalue, res_df$padj, n_case, n_ctrl, accession,
                  "single_cell", detect_id_space_r(rownames(mat)), "log2FC")
}
