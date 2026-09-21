library(testthat)

# single_cell_da.R cannot be source()d here: its top-level block loads DESeq2,
# which is absent under the R the CI job installs. These are source-level
# regression guards over the conventions the arm shares with the bulk arm --
# the only things a wrong edit could silently break, because a pseudobulk
# effect table is the right shape whatever standard error it carries.

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
SRC <- readLines(file.path(REPO_ROOT, "R", "single_cell_da.R"), warn = FALSE)
CODE <- SRC[!grepl("^\\s*#", SRC)]

STEP05 <- readLines(file.path(REPO_ROOT, "pipeline", "05_per_study_da.R"),
                    warn = FALSE)
STEP06 <- readLines(file.path(REPO_ROOT, "pipeline", "06_meta_analysis.R"),
                    warn = FALSE)


test_that("apeglm shrinkage never returns to the single-cell arm either", {
  # The same defect as in the bulk arm: a posterior SD of ~0.001 against a true
  # SE of 0.10-0.25 would collapse the random-effects weights. It matters more
  # here, since the point of the arm is to be pooled against the bulk estimates.
  expect_false(any(grepl("lfcShrink", CODE)))
  expect_false(any(grepl("apeglm", CODE)))
})

test_that("the SE comes from the unshrunken results() object", {
  expect_true(any(grepl("results\\(dds", CODE)))
  expect_true(any(grepl("res_df\\$lfcSE", CODE)))
})

test_that("the effect table is built by the shared constructor", {
  # A hand-rolled data.frame here would drift from the bulk arm's columns, and
  # step 06 would pool the two on a schema that only looks the same.
  expect_true(any(grepl("build_effect_df\\(", CODE)))
  expect_true(any(grepl("da_utils\\.R", SRC)))
})

test_that("the arm is labelled single_cell, not transcriptomics", {
  # Mislabelling the modality column is how the arm's rows would end up
  # indistinguishable from bulk rows in any table that concatenates them.
  expect_true(any(grepl('"single_cell"', CODE)))
  expect_false(any(grepl('"transcriptomics"', CODE)))
})

test_that("the arm counts are taken per subject, not per GSM", {
  # Two libraries of one animal share a subject id. Reading the assignment per
  # GSM would enter that animal twice and inflate n by the number of repeats.
  expect_true(any(grepl("by_subject", CODE)))
  expect_true(any(grepl("groups\\[, c\\(\"subject\", \"group\"\\)\\]", CODE)))
})

test_that("the minimum samples per group is honoured", {
  expect_true(any(grepl("min_samples_per_group", CODE)))
})

test_that("groups are read from data/interim, never from the shared geo_cache", {
  # Every accession in this arm is also in the bulk arm's screening sheet. A
  # groups file in data/raw/geo_cache/ is what the bulk transcriptomics DA
  # reads, so writing one there would wire the two arms together.
  expect_true(any(grepl("single_cell", CODE)))
  expect_false(any(grepl("geo_cache", CODE)))
})

test_that("nothing in the arm downloads from GEO at analysis time", {
  # Step 19 resolved four deposit formats already; a getGEO here would go
  # around it and read whatever the series matrix happens to carry.
  expect_false(any(grepl("getGEO", CODE)))
})


test_that("step 05 dispatches single_cell and keeps it out of ALL_MODALITIES", {
  all_line <- grep("^ALL_MODALITIES <-", STEP05, value = TRUE)
  expect_length(all_line, 1)
  expect_false(grepl("single_cell", all_line))
  expect_true(any(grepl("single_cell = run_single_cell_da\\(", STEP05)))
  expect_true(any(grepl("single_cell_da\\.R", STEP05)))
})

test_that("step 06 keeps single_cell out of ALL_MODALITIES", {
  all_line <- grep("^ALL_MODALITIES <-", STEP06, value = TRUE)
  expect_length(all_line, 1)
  expect_false(grepl("single_cell", all_line))
})

test_that("both steps refuse a modality they cannot dispatch", {
  # Without this a typo produces an empty results directory and a green run.
  for (src in list(STEP05, STEP06)) {
    expect_true(any(grepl("unknown modality", src)))
    expect_true(any(grepl("EXTRA_MODALITIES", src)))
  }
})
