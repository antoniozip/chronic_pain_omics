library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "da_utils.R"))
source(file.path(REPO_ROOT, "R", "genomics_da.R"))


test_that("parse_sample_count extracts integer from descriptor", {
  expect_equal(parse_sample_count("4,326 European ancestry"), 4326L)
  expect_equal(parse_sample_count("1000 individuals"), 1000L)
  expect_equal(parse_sample_count(""), 0L)
  expect_equal(parse_sample_count(NULL), 0L)
})

test_that("parse_sample_count handles no digits gracefully", {
  expect_equal(parse_sample_count("N/A"), 0L)
  expect_equal(parse_sample_count("unknown"), 0L)
})

test_that("standardize_gwas drops rows with missing rsid", {
  df <- data.frame(
    rsid        = c("rs123", NA, "rs456"),
    beta_or     = c(0.05, 0.1, -0.03),
    or_per_copy = c(1.2, 1.1, 0.9),
    pval        = c(5e-8, 0.01, 1e-6),
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "beta")
  expect_equal(nrow(result), 2)
  expect_true(all(!is.na(result$feature_id)))
})

test_that("standardize_gwas uses log(OR) for OR effect unit", {
  df <- data.frame(
    rsid = "rs123", beta_or = 0.1, or_per_copy = 1.5, pval = 5e-8,
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "OR")
  expect_equal(result$effect_size, log(1.5), tolerance = 1e-6)
})

test_that("standardize_gwas uses beta for beta effect unit", {
  df <- data.frame(
    rsid = "rs123", beta_or = 0.42, or_per_copy = NA_real_, pval = 5e-8,
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "beta")
  expect_equal(result$effect_size, 0.42, tolerance = 1e-6)
})
