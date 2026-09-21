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

test_that("isoform rows are deduplicated keeping the minimum p-value", {
  rows <- base_rows()
  rows <- rbind(rows, data.frame(
    gene = "Bdnf", `log2(fold_change)` = 3.0, p_value = 0.001,
    q_value = 0.002, test_stat = 9.9, status = "OK",
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  df <- run_cuffdiff_output(write_cuffdiff(rows), "GSE1", list())

  expect_equal(sum(df$feature_id == "Bdnf"), 1)
  expect_equal(df$pval[df$feature_id == "Bdnf"], 0.001)
  expect_equal(df$effect_size[df$feature_id == "Bdnf"], 3.0)
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
