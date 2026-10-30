library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "da_utils.R"))
source(file.path(REPO_ROOT, "R", "metabolomics_da.R"))


# ---------------------------------------------------------------------------
# compute_hedges_g
# ---------------------------------------------------------------------------

test_that("compute_hedges_g returns numeric vector with correct length", {
  set.seed(42)
  mat   <- matrix(rnorm(40), nrow = 4, ncol = 10,
                  dimnames = list(paste0("M", 1:4), paste0("S", 1:10)))
  group <- factor(c(rep("case", 5), rep("control", 5)))
  result <- compute_hedges_g(mat, group)
  expect_length(result, 4)
  expect_type(result, "double")
})

test_that("compute_hedges_g is near zero for same-distribution groups", {
  set.seed(1)
  mat   <- matrix(rnorm(200), nrow = 10, ncol = 20,
                  dimnames = list(paste0("M", 1:10), paste0("S", 1:20)))
  group <- factor(rep(c("case", "control"), each = 10))
  g     <- compute_hedges_g(mat, group)
  # With equal distributions, median |g| should be small
  expect_lt(median(abs(g), na.rm = TRUE), 1.0)
})

test_that("compute_hedges_g is large for well-separated groups", {
  set.seed(7)
  case_vals <- matrix(rnorm(50, mean = 5, sd = 0.5), nrow = 5)
  ctrl_vals <- matrix(rnorm(50, mean = 0, sd = 0.5), nrow = 5)
  mat   <- cbind(case_vals, ctrl_vals)
  dimnames(mat) <- list(paste0("M", 1:5), paste0("S", 1:20))
  group <- factor(c(rep("case", 10), rep("control", 10)))
  g     <- compute_hedges_g(mat, group)
  expect_true(all(abs(g) > 3, na.rm = TRUE))
})

test_that("compute_hedges_g handles NA values", {
  mat        <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8), nrow = 2)
  dimnames(mat) <- list(c("M1", "M2"), paste0("S", 1:4))
  group <- factor(c("case", "case", "control", "control"))
  result <- compute_hedges_g(mat, group)
  expect_length(result, 2)
  # Should not throw even with NAs
})


# ---------------------------------------------------------------------------
# the SE must belong to the effect that is reported
#
# The backend reports Hedges' g but carried the limma SE of log2FC. 06/06b/06g
# weight on 1/se^2 and pool effect_size, so every metabolomic feature was
# weighted by the precision of a different quantity on a different scale.
# ---------------------------------------------------------------------------

test_that("compute_hedges_g returns the standard error of g on request", {
  set.seed(11)
  mat   <- matrix(rnorm(60, mean = rep(c(2, 0), each = 30)), nrow = 3, ncol = 20)
  group <- factor(rep(c("case", "control"), each = 10))
  n1 <- sum(group == "case"); n0 <- sum(group == "control")

  out <- compute_hedges_g(mat, group, return_se = TRUE)
  expect_type(out, "list")
  expect_named(out, c("g", "se"))
  expect_equal(out$g, compute_hedges_g(mat, group))

  expected <- sqrt((n1 + n0) / (n1 * n0) + out$g^2 / (2 * (n1 + n0)))
  expect_equal(unname(out$se), unname(expected))
})

test_that("the default return of compute_hedges_g is still the g vector", {
  mat   <- matrix(rnorm(20), nrow = 2, ncol = 10)
  group <- factor(rep(c("case", "control"), each = 5))
  expect_type(compute_hedges_g(mat, group), "double")
})

test_that("compute_hedges_g's SE is not the SE of log2FC", {
  set.seed(3)
  mat   <- matrix(rnorm(40, mean = 1), nrow = 4, ncol = 10)
  group <- factor(rep(c("case", "control"), each = 5))
  out   <- compute_hedges_g(mat, group, return_se = TRUE)
  # limma's SE is on the log2 concentration scale and is much smaller for
  # well-measured features; the two must not be interchangeable.
  expect_true(all(out$se > 0.3))
})


# ---------------------------------------------------------------------------
# run_metabolomics_da (smoke test with synthetic data)
# ---------------------------------------------------------------------------

test_that("run_metabolomics_da returns NULL for insufficient groups", {
  set.seed(42)
  mat <- matrix(rnorm(40), nrow = 4, ncol = 10,
                dimnames = list(paste0("M", 1:4), paste0("S", 1:10)))
  # Only case samples — no controls
  meta <- data.frame(group = rep("case", 10), row.names = paste0("S", 1:10))
  cfg  <- list(da = list(min_samples_per_group = 3,
                         case_labels = "case", control_labels = "control"))
  result <- run_metabolomics_da("MTBLS_TEST", mat, meta, "ChEBI", "metabolomics", cfg)
  expect_null(result)
})

test_that("run_metabolomics_da returns effect df with correct columns", {
  set.seed(42)
  n <- 12
  mat <- matrix(c(rnorm(6, 5), rnorm(6, 0)), nrow = 3, ncol = n,
                dimnames = list(c("M1", "M2", "M3"), paste0("S", 1:n)))
  meta <- data.frame(
    group = rep(c("case", "control"), each = 6),
    row.names = paste0("S", 1:n),
    stringsAsFactors = FALSE
  )
  cfg <- list(da = list(
    min_samples_per_group = 3,
    case_labels  = "case",
    control_labels = "control"
  ))
  result <- run_metabolomics_da("MTBLS_TEST", mat, meta, "ChEBI", "metabolomics", cfg)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("feature_id", "effect_size", "se", "pval") %in% colnames(result)))
  expect_equal(nrow(result), 3)
  expect_equal(result$id_space[1], "ChEBI")
  expect_equal(result$effect_unit[1], "SMD")
  # SMD is only truthful if the SE is the SE of g. It was limma's abs(logFC/t)
  # — the precision of a quantity on a different scale.
  n1 <- sum(meta$group == "case")
  n0 <- sum(meta$group == "control")
  expected_se <- sqrt((n1 + n0) / (n1 * n0) +
                        result$effect_size^2 / (2 * (n1 + n0)))
  expect_equal(result$se, unname(expected_se), tolerance = 1e-8)
})
