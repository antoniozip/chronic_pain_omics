library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "se_utils.R"))


test_that("se_from_t computes |logFC / t|", {
  expect_equal(se_from_t(2, 4), 0.5)
  expect_equal(se_from_t(-1.5, 3), 0.5)          # sign-independent
  expect_equal(se_from_t(c(2, -1.5), c(4, 3)), c(0.5, 0.5))
})

test_that("se_from_t is not the qt-based formula", {
  # Guard against reverting to logFC / qt(0.975, df), which inflates SE.
  expect_false(isTRUE(all.equal(se_from_t(2, 4), 2 / qt(0.975, 10))))
})

test_that("se_from_pvalue inverts the two-sided Normal approximation", {
  expect_equal(se_from_pvalue(1, 0.05), 1 / qnorm(0.975), tolerance = 1e-9)
  expect_equal(round(se_from_pvalue(1, 0.05), 7), 0.5102135)
  expect_equal(round(se_from_pvalue(2, 0.01), 7), 0.7764490)
})

test_that("se_from_pvalue clamps p to avoid infinite z", {
  s <- se_from_pvalue(1, 0)
  expect_true(is.finite(s))
  expect_equal(round(s, 7), 0.1546324)
})

test_that("se_from_pvalue returns NA for degenerate inputs", {
  expect_true(is.na(se_from_pvalue(0, 0.05)))   # zero effect
  expect_true(is.na(se_from_pvalue(1, 1)))      # p=1 => z=0
  expect_true(is.na(se_from_pvalue(NA_real_, 0.05)))
})

test_that("needs_log2 uses a Q75 threshold of 50, not 100", {
  # GSE10238 has Q75 = 92.2 and MUST be transformed.
  expect_true(needs_log2(c(rep(1, 3), 92.2, 200)))
  expect_false(needs_log2(c(1, 2, 3, 4, 5)))
})

test_that("maybe_log2 transforms only when needed and floors at 1", {
  linear <- c(1, 4, 256, 1024)
  expect_equal(maybe_log2(linear), log2(pmax(linear, 1)))

  already_log <- c(1, 2, 3, 4)
  expect_equal(maybe_log2(already_log), already_log)

  expect_equal(maybe_log2(c(0, 4, 256, 1024))[1], 0)  # log2(pmax(0,1)) == 0
})
