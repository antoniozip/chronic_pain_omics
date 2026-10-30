library(testthat)

# 06d_sni_tissue_stratified.R is a runnable script: sourcing it reads results/
# and fits the whole SNI moderator set. Only its pure helper definitions are
# evaluated here, the same way test-normalize-ids.R lifts individual functions.

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))

SNI_SCRIPT <- file.path(REPO_ROOT, "pipeline", "06d_sni_tissue_stratified.R")

for (expr in parse(SNI_SCRIPT)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, envir = environment())
  }
}


# ---------------------------------------------------------------------------
# Y-09: the BH family must be pre-specified
#
# The tested set used to be `head(sni_pooled_idx$feature_id, 200)` after
# `order(-abs(yi))`: the family was chosen on effect magnitude, so the 200 kept
# features were enriched for small p-values and the BH-adjusted count was
# anti-conservative relative to correcting over every candidate feature. The
# family is now every feature that qualifies on the data-independent rule
# (>= 2 tissues and >= 3 studies).
# ---------------------------------------------------------------------------

test_that("the tissue-moderator family is every qualifying feature, not the top |yi|", {
  df <- data.frame(
    feature_id  = c(rep("A", 3), rep("B", 3), rep("C", 2)),
    study_id    = c("s1", "s2", "s3", "s1", "s2", "s3", "s1", "s2"),
    tissue      = c("DRG", "Spinal_cord", "mPFC",
                    "DRG", "Spinal_cord", "DRG",
                    "DRG", "DRG"),
    effect_size = c(10, 10, 10, -0.1, -0.1, -0.1, 0, 0),
    se          = 0.1,
    stringsAsFactors = FALSE
  )

  # A: 3 tissues / 3 studies. B: 2 tissues / 3 studies. C: 1 tissue.
  expect_setequal(sni_tissue_family(df), c("A", "B"))

  # Membership must not move when the effects do: the rule is on the design of
  # the stratum, not on the result.
  df$effect_size <- rev(df$effect_size)
  expect_setequal(sni_tissue_family(df), c("A", "B"))
})

test_that("the family is not capped at the 200 largest effects", {
  n <- 250
  df <- data.frame(
    feature_id  = rep(sprintf("G%03d", seq_len(n)), each = 3),
    study_id    = rep(c("s1", "s2", "s3"), n),
    tissue      = rep(c("DRG", "Spinal_cord", "mPFC"), n),
    effect_size = rnorm(3 * n),
    se          = 0.2,
    stringsAsFactors = FALSE
  )
  expect_length(sni_tissue_family(df), n)
})


# ---------------------------------------------------------------------------
# Y-08: Q_between_df is the moderator df, not the coefficient count
#
# `fit$m` counts the moderator model's coefficients, intercept included, so a
# three-tissue test reported df 3 where the QM test has 2. R/meta_utils.R
# moderator_df() exists for exactly this and 06b already uses it; 06d wrote
# fit$m into sni_tissue_meta.csv, a published supplementary input.
# ---------------------------------------------------------------------------

test_that("the moderator df is the QM df (g-1), not the coefficient count", {
  sub <- data.frame(
    effect_size = c(0.5, 0.6, 0.1, 0.2, 0.8, 0.7),
    se          = rep(0.1, 6),
    tissue      = rep(c("DRG", "Spinal_cord", "mPFC"), each = 2),
    stringsAsFactors = FALSE
  )

  row <- sni_tissue_moderator(sub, "Gene1")
  fit <- metafor::rma(yi = effect_size, sei = se, mods = ~ tissue,
                      data = sub, method = "REML")

  expect_equal(row$Q_between_df, moderator_df(fit))
  expect_equal(row$Q_between_df, 2)      # three tissues -> g - 1
  expect_equal(fit$p, 3)                 # coefficients, intercept included
  expect_false(isTRUE(all.equal(row$Q_between_df, fit$p)))
})

test_that("a feature with fewer than two tissues has no moderator test", {
  sub <- data.frame(
    effect_size = c(0.5, 0.6, 0.1),
    se = rep(0.1, 3),
    tissue = rep("DRG", 3),
    stringsAsFactors = FALSE
  )
  expect_null(sni_tissue_moderator(sub, "Gene1"))
})

# The gate used to be `fit$m < 2`. Under metafor's actual definition `$m` is
# the number of moderator coefficients, so a two-tissue feature (one moderator,
# df 1) was silently not tested at all -- while the family rule admits every
# feature seen in >= 2 tissues. The gate is now "at least one moderator df".
test_that("a two-tissue feature is tested with one moderator df", {
  sub <- data.frame(
    effect_size = c(0.5, 0.6, 0.1, 0.2, 0.9, 0.8),
    se = rep(0.1, 6),
    tissue = rep(c("DRG", "Spinal_cord"), each = 3),
    stringsAsFactors = FALSE
  )
  row <- sni_tissue_moderator(sub, "Gene1")
  expect_equal(row$Q_between_df, 1)
  expect_equal(row$n_tissues, 2)
})
