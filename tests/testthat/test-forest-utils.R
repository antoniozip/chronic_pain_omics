library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))
source(file.path(REPO_ROOT, "R", "forest_utils.R"))


#' Page dimensions of a PDF, in inches, read from its MediaBox entries.
page_dimensions <- function(path) {
  raw <- readBin(path, "raw", file.size(path))
  # Matched on bytes: a PDF is not valid UTF-8, and a regex over it otherwise
  # fails on the first stream that happens not to decode.
  txt <- rawToChar(raw[raw != as.raw(0)])
  boxes <- regmatches(txt, gregexpr("/MediaBox *\\[[^]]*\\]", txt,
                                    useBytes = TRUE))[[1]]
  dims <- lapply(boxes, function(b) {
    nums <- as.numeric(strsplit(gsub("[^0-9. ]", " ", b), " +")[[1]])
    nums <- nums[!is.na(nums)]
    c(width = nums[3], height = nums[4]) / 72
  })
  data.frame(width = vapply(dims, `[[`, numeric(1), "width"),
             height = vapply(dims, `[[`, numeric(1), "height"))
}

#' Page heights of a PDF, in inches, read from its MediaBox entries.
page_heights <- function(path) {
  page_dimensions(path)$height
}

#' Per-study rows for `n` studies of one feature, with usable standard errors.
feature_rows <- function(feature_id, n) {
  data.frame(
    feature_id  = feature_id,
    study_id    = paste0("GSE", seq_len(n)),
    effect_size = seq(-0.5, 0.5, length.out = n),
    se          = rep(0.2, n),
    effect_unit = "log2FC",
    stringsAsFactors = FALSE
  )
}


test_that("forest_page_height grows with k above a floor", {
  expect_equal(forest_page_height(1), 3.2)
  expect_equal(forest_page_height(4), 3.4)
  expect_equal(forest_page_height(90), 33.5)
  expect_gt(forest_page_height(50), forest_page_height(20))
})


test_that("one feature gets a page sized to its own panel", {
  out <- file.path(tempdir(), "one.pdf")
  on.exit(unlink(out), add = TRUE)
  write_forest_pdf(out, feature_rows("A", 5), "A")

  expect_true(file.exists(out))
  expect_equal(page_heights(out), forest_page_height(5), tolerance = 0.01)
})


test_that("a set sizes every page to its own panel, not to the deepest", {
  # The failure this guards: one pdf() device fixes one page size, so the
  # whole set took the deepest panel's page and forest() stretched a
  # five-study panel down 34 inches.
  skip_if(!nzchar(Sys.which("pdfunite")) && !nzchar(Sys.which("gs")),
          "no pdfunite or gs to join the pages")
  effects <- rbind(feature_rows("shallow", 3), feature_rows("deep", 40))
  out <- file.path(tempdir(), "set.pdf")
  on.exit(unlink(out), add = TRUE)
  write_forest_pdf(out, effects, c("shallow", "deep"))

  heights <- page_heights(out)
  expect_length(heights, 2)
  expect_equal(heights, c(forest_page_height(3), forest_page_height(40)),
               tolerance = 0.01)
})


test_that("a feature with too few usable rows contributes no page", {
  skip_if(!nzchar(Sys.which("pdfunite")) && !nzchar(Sys.which("gs")),
          "no pdfunite or gs to join the pages")
  effects <- rbind(feature_rows("drawable", 4), feature_rows("single", 1))
  out <- file.path(tempdir(), "partial.pdf")
  on.exit(unlink(out), add = TRUE)
  write_forest_pdf(out, effects, c("drawable", "single"))

  expect_equal(page_heights(out), forest_page_height(4), tolerance = 0.01)
})


test_that("a feature absent from the effects table is skipped", {
  out <- file.path(tempdir(), "absent.pdf")
  on.exit(unlink(out), add = TRUE)
  expect_null(write_forest_pdf(out, feature_rows("A", 4), "B"))
  expect_false(file.exists(out))
})


# Figure 2 grew from one feature to the ten most significant, and the ten are
# shallow (3-7 studies each), so they fit on one page as a grid. The
# supplementary set keeps one page per feature, where the deepest panel can be
# 90 studies and no grid would be readable.
test_that("a forest grid draws every feature on one page", {
  out <- file.path(tempdir(), "grid.pdf")
  on.exit(unlink(out), add = TRUE)
  effects <- rbind(feature_rows("A", 3), feature_rows("B", 5),
                   feature_rows("C", 4), feature_rows("D", 7))
  write_forest_grid_pdf(out, effects, c("A", "B", "C", "D"),
                        ncol = 2, width = 10, panel_height = 2.8)

  dims <- page_dimensions(out)
  expect_equal(nrow(dims), 1)
  expect_equal(dims$width, 10, tolerance = 0.01)
  expect_equal(dims$height, 2 * 2.8, tolerance = 0.01)
})


test_that("a grid sized its page from the features it can draw", {
  out <- file.path(tempdir(), "grid-partial.pdf")
  on.exit(unlink(out), add = TRUE)
  # Three features requested, one absent and one single-study: only two panels
  # are drawable, which is one row at ncol = 2 rather than the requested two.
  effects <- rbind(feature_rows("A", 3), feature_rows("B", 4),
                   feature_rows("single", 1))
  write_forest_grid_pdf(out, effects, c("A", "B", "single"), ncol = 2,
                        width = 10, panel_height = 2.8)

  expect_equal(page_heights(out), 2.8, tolerance = 0.01)
})


test_that("a grid with no drawable feature writes no file", {
  out <- file.path(tempdir(), "grid-empty.pdf")
  on.exit(unlink(out), add = TRUE)
  expect_null(write_forest_grid_pdf(out, feature_rows("A", 1), c("A", "B")))
  expect_false(file.exists(out))
})
