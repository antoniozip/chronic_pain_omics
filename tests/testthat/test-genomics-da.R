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


# ---------------------------------------------------------------------------
# se and padj
#
# Every row used to be emitted with se = NA_real_ and a comment about an
# imputation from the CI that was never written. 06/06b/06g all drop
# `is.na(se) | se <= 0`, so a genomics run produced a *successful* empty pool.
# And padj was the raw p-value, so the effect table claimed padj == pval.
# ---------------------------------------------------------------------------

test_that("standardize_gwas emits a real se derived from the effect and p", {
  df <- data.frame(
    rsid = "rs123", beta_or = 0.05, or_per_copy = NA_real_, pval = 5e-8,
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "beta")
  expect_equal(result$se, 0.05 / qnorm(1 - 5e-8 / 2), tolerance = 1e-10)
})

test_that("standardize_gwas derives se for the log(OR) effect as well", {
  df <- data.frame(
    rsid = "rs123", beta_or = NA_real_, or_per_copy = 1.5, pval = 1e-6,
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "OR")
  expect_equal(result$se, abs(log(1.5)) / qnorm(1 - 1e-6 / 2), tolerance = 1e-10)
})

test_that("standardize_gwas drops rows that cannot yield a usable se", {
  df <- data.frame(
    rsid = c("rs_ok", "rs_p1", "rs_nop", "rs_zero"),
    beta_or = c(0.2, 0.2, 0.2, 0.0),
    or_per_copy = NA_real_,
    pval = c(1e-4, 1.0, NA_real_, 1e-4),
    stringsAsFactors = FALSE
  )
  result <- standardize_gwas(df, "beta")
  # p == 1 gives se = 0, p = NA gives se = NA, a zero effect gives se = 0:
  # all three are dropped rather than emitted with an unusable weight.
  expect_equal(result$feature_id, "rs_ok")
  expect_true(all(is.finite(result$se) & result$se > 0))
})

test_that("the per-study padj is BH through the shared FDR helper, never raw p", {
  df <- data.frame(pval = c(0.001, 0.01, 0.02, 0.5))
  out <- gwas_padj(df)
  expect_equal(out$padj, p.adjust(df$pval, method = "BH"))
  expect_false(isTRUE(all.equal(out$padj, out$pval)))
})
