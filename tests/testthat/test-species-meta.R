library(testthat)

# 06g_species_meta.R cannot be source()d here: its top-level block loads
# metafor and reads the per-study tree. As elsewhere in this suite, lift the
# individual functions out of the source.

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
SRC <- readLines(file.path(REPO_ROOT, "pipeline", "06g_species_meta.R"), warn = FALSE)
CODE <- SRC[!grepl("^\\s*#", SRC)]

lift <- function(name) {
  at <- grep(paste0("^", name, " <- function"), SRC)
  stopifnot(length(at) == 1)
  ends <- grep("^\\}", SRC)
  eval(parse(text = paste(SRC[at:min(ends[ends > at])], collapse = "\n")),
       envir = parent.frame())
}


# ---------------------------------------------------------------------------
# --exclude-studies
#
# The 2026-08-29 cross-species reversal needs a tissue-matched sensitivity
# run: the human pool grew from 6 studies to 24, 14 of them endometrial, and
# concordance moved from significantly below chance to significantly above it.
# Answering "is this composition?" means dropping those studies -- which must
# not be done through superseded_studies.csv, because that registry records a
# claim about sample independence and would be left asserting something false.
# ---------------------------------------------------------------------------

test_that("read_excluded_studies returns nothing when no file is given", {
  lift("read_excluded_studies")
  expect_equal(read_excluded_studies(NA_character_), character(0))
  expect_equal(read_excluded_studies(""), character(0))
})

test_that("read_excluded_studies returns nothing for a path that does not exist", {
  lift("read_excluded_studies")
  expect_equal(read_excluded_studies("/nonexistent/exclusions.csv"), character(0))
})

test_that("read_excluded_studies reads the study_id column", {
  lift("read_excluded_studies")
  f <- withr::local_tempfile(fileext = ".csv")
  write.csv(data.frame(study_id = c("GSE1", "GSE2"), reason = c("a", "b")),
            f, row.names = FALSE)
  expect_equal(read_excluded_studies(f), c("GSE1", "GSE2"))
})

test_that("read_excluded_studies refuses a file without a study_id column", {
  # A silently ignored exclusion file would produce a sensitivity result
  # identical to the main one, which reads as "composition makes no
  # difference" -- the worst possible failure for this analysis.
  lift("read_excluded_studies")
  f <- withr::local_tempfile(fileext = ".csv")
  write.csv(data.frame(accession = "GSE1"), f, row.names = FALSE)
  expect_error(read_excluded_studies(f), "study_id")
})

test_that("the exclusion is applied where studies are selected, not after pooling", {
  at_filter <- grep("EXCLUDED_STUDIES", CODE)
  at_pool <- grep("rma\\(", CODE)
  expect_gt(length(at_filter), 0)
  expect_lt(min(at_filter[at_filter > 100]), min(at_pool))
})

test_that("the output directory can be redirected away from the published pools", {
  expect_true(any(grepl("OUT_DIR_ARG", CODE)))
  expect_true(any(grepl('match\\("--out-dir", args\\)', CODE)))
})
