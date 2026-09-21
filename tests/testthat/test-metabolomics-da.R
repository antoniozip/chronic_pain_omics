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
})
