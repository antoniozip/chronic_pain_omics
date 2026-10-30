library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "da_utils.R"))
source(file.path(REPO_ROOT, "R", "se_utils.R"))
source(file.path(REPO_ROOT, "R", "cuffdiff_da.R"))


write_cuffdiff <- function(rows) {
  path <- tempfile(fileext = ".diff")
  write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

base_rows <- function() {
  data.frame(
    gene                = c("Bdnf", "Cd36"),
    `log2(fold_change)` = c(1.0, 2.0),
    p_value             = c(0.05, 0.01),
    q_value             = c(0.10, 0.02),
    test_stat           = c(3.3, 7.7),   # deliberately NOT logFC/SE
    status              = c("OK", "OK"),
    sample_1            = c("control", "control"),
    sample_2            = c("case", "case"),
    check.names         = FALSE,
    stringsAsFactors    = FALSE
  )
}


test_that("SE is derived from the p-value, not from test_stat", {
  df <- run_cuffdiff_output(write_cuffdiff(base_rows()), "GSE1", list())
  se <- df$se[df$feature_id == "Bdnf"]

  expect_equal(round(se, 7), 0.5102135)          # 1.0 / qnorm(0.975)
  expect_false(isTRUE(all.equal(se, 1.0 / 3.3))) # would be the test_stat bug
})

test_that("rows with unusable status are dropped", {
  rows <- base_rows()
  rows$status <- c("OK", "FAIL")
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(nrow(df), 1)
  expect_equal(df$feature_id, "Bdnf")
})

test_that("placeholder gene names are dropped", {
  rows <- base_rows()
  rows$gene <- c("-", "Cd36")
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(df$feature_id, "Cd36")
})

test_that("isoform rows are deduplicated to one row per gene", {
  rows <- base_rows()
  rows <- rbind(rows, data.frame(
    gene = "Bdnf", `log2(fold_change)` = 3.0, p_value = 0.001,
    q_value = 0.002, test_stat = 9.9, status = "OK",
    sample_1 = "control", sample_2 = "case",
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())

  expect_equal(sum(df$feature_id == "Bdnf"), 1)
})

# The rule used to be `order(pval)` + `duplicated(gene)`: the retained isoform
# was the most significant one, which is selection on the outcome (05b removed
# the identical rule from aggregate_to_gene for the same reason -- the minimum
# of n p-values is not uniform under the null and the effect is biased away
# from zero). The rule is now largest |log2FC|, pre-specified and independent
# of the p-value.
test_that("isoforms are not chosen by smallest p-value", {
  rows <- base_rows()
  rows[1, "log2(fold_change)"] <- 2.5      # Bdnf: strong effect, p = 0.001
  rows[1, "p_value"]           <- 0.001
  rows <- rbind(rows, data.frame(
    gene = "Bdnf", `log2(fold_change)` = 3.0,   # larger effect, weaker p
    p_value = 0.05, q_value = 0.10, test_stat = 1.9, status = "OK",
    sample_1 = "control", sample_2 = "case",
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())

  bdnf <- df[df$feature_id == "Bdnf", ]
  expect_equal(nrow(bdnf), 1)
  expect_equal(bdnf$effect_size, 3.0)
  expect_equal(bdnf$pval, 0.05)            # not the min-p row
})

# ---------------------------------------------------------------------------
# Y-04: Cuffdiff's log2(fold_change) is log2(sample_2 / sample_1), so the sign
# is established by the file's own sample names. Every other backend fixes
# case-as-positive (limma coef "groupcase", m1 - m0, contrast case-control);
# a flipped Cuffdiff study cancels the others on a shared feature instead of
# reinforcing them, and se_from_pvalue() is symmetric so nothing notices.
# ---------------------------------------------------------------------------

test_that("a control-as-sample_1 file keeps its sign (case minus control)", {
  rows <- base_rows()   # sample_1 = control, sample_2 = case
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(df$effect_size[df$feature_id == "Bdnf"], 1.0)
})

test_that("a case-as-sample_1 file is flipped to case minus control", {
  rows <- base_rows()
  rows$sample_1 <- "case"
  rows$sample_2 <- "control"
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(df$effect_size[df$feature_id == "Bdnf"], -1.0)
  expect_equal(df$effect_size[df$feature_id == "Cd36"], -2.0)
})

test_that("the contrast is resolved from the configured arm vocabulary too", {
  rows <- base_rows()
  rows$sample_1 <- "healthy control"     # control arm is the denominator
  rows$sample_2 <- "fibromyalgia"        # case arm is the numerator: keep sign
  cfg <- list(da = list(case_labels = "fibromyalgia",
                        control_labels = c("healthy", "control")))
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", cfg)
  expect_equal(df$effect_size[df$feature_id == "Bdnf"], 1.0)
})

test_that("a file whose sample names match neither arm is refused", {
  rows <- base_rows()
  rows$sample_1 <- "treated"
  rows$sample_2 <- "vehicle"
  expect_error(run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list()),
               "contrast direction")
})

test_that("a file without sample_1/sample_2 cannot establish the sign", {
  rows <- base_rows()
  rows$sample_1 <- NULL
  rows$sample_2 <- NULL
  expect_error(run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list()),
               "sample_1")
})

test_that("zero-effect rows are dropped because SE is undefined", {
  rows <- base_rows()
  rows[["log2(fold_change)"]] <- c(0.0, 2.0)
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(df$feature_id, "Cd36")
})

test_that("missing required columns yield NULL", {
  rows <- base_rows()
  rows$test_stat <- NULL
  expect_null(run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list()))
})

test_that("output carries the build_effect_df contract", {
  df <- run_cuffdiff_output(write_cuffdiff(base_rows()), "GSE42", list())
  expect_true(all(c("feature_id", "effect_size", "se", "pval", "padj",
                    "study_id", "modality", "id_space", "effect_unit")
                  %in% colnames(df)))
  expect_equal(unique(df$study_id), "GSE42")
  expect_equal(unique(df$modality), "transcriptomics")
  expect_equal(unique(df$effect_unit), "log2FC")
  expect_type(df$feature_id, "character")
})

test_that("NULL is returned when every row is filtered out", {
  rows <- base_rows()
  rows[["log2(fold_change)"]] <- c(0, 0)   # se_from_pvalue returns NA for zero effect
  expect_null(run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list()))
})

test_that("NOTEST status rows are retained", {
  rows <- base_rows()
  rows$status <- c("OK", "NOTEST")
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(nrow(df), 2)
})

test_that("empty string gene names are dropped", {
  rows <- base_rows()
  rows$gene <- c("", "Cd36")
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())
  expect_equal(df$feature_id, "Cd36")
})

test_that("nonexistent file path returns NULL", {
  expect_null(suppressWarnings(suppressMessages(run_cuffdiff_output("/nonexistent/path/to/file.diff", "GSE1", list()))))
})
