library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "symbol_utils.R"))


# ---------------------------------------------------------------------------
# symbol_case_style
# ---------------------------------------------------------------------------

test_that("case style follows the majority convention in the table", {
  expect_equal(symbol_case_style(c("TNF", "ACTB", "GAPDH", "IL6")), "upper")
  expect_equal(symbol_case_style(c("Tnf", "Actb", "Gapdh", "Il6")), "title")
})

test_that("numeric probe ids do not decide the case style", {
  # Rodent tables carry many bare numeric ids; they carry no case information
  # and must not outvote the gene symbols.
  expect_equal(symbol_case_style(c("10788858", "10930610", "Tnf", "Actb")), "title")
})


# ---------------------------------------------------------------------------
# Excel date corruption
# ---------------------------------------------------------------------------

test_that("MARCH family dates are repaired to the MARCHF symbol", {
  expect_equal(canonicalize_gene_symbols("1-Mar", "title"), "Marchf1")
  expect_equal(canonicalize_gene_symbols("11-Mar", "title"), "Marchf11")
  expect_equal(canonicalize_gene_symbols("1-Mar", "upper"), "MARCHF1")
})

test_that("septin dates are repaired to the SEPTIN symbol", {
  expect_equal(canonicalize_gene_symbols("1-Sep", "title"), "Septin1")
  expect_equal(canonicalize_gene_symbols("14-Sep", "title"), "Septin14")
  expect_equal(canonicalize_gene_symbols("2-Sep", "upper"), "SEPTIN2")
})

test_that("15-Sep is SEP15/SELENOF, not a septin", {
  # There is no SEPT15. 15-Sep is the Excel rendering of SEP15, the 15 kDa
  # selenoprotein, now SELENOF. Mapping it to Septin15 would invent a gene.
  expect_equal(canonicalize_gene_symbols("15-Sep", "title"), "Selenof")
  expect_equal(canonicalize_gene_symbols("15-Sep", "upper"), "SELENOF")
})

test_that("the month-first spelling is repaired too", {
  expect_equal(canonicalize_gene_symbols("Mar-3", "title"), "Marchf3")
  expect_equal(canonicalize_gene_symbols("Sep-9", "title"), "Septin9")
})

test_that("dates outside the known families are left alone", {
  # SEPT13 does not exist, so 13-Sep cannot be attributed to a gene.
  expect_equal(canonicalize_gene_symbols("13-Sep", "title"), "13-Sep")
  expect_equal(canonicalize_gene_symbols("12-Mar", "title"), "12-Mar")
  expect_equal(canonicalize_gene_symbols("5-Dec", "title"), "5-Dec")
})


# ---------------------------------------------------------------------------
# Legacy aliases and casing
# ---------------------------------------------------------------------------

test_that("legacy MARCH and SEPT names map to current symbols", {
  expect_equal(canonicalize_gene_symbols("March1", "title"), "Marchf1")
  expect_equal(canonicalize_gene_symbols("MARCH1", "upper"), "MARCHF1")
  expect_equal(canonicalize_gene_symbols("Sept1", "title"), "Septin1")
  expect_equal(canonicalize_gene_symbols("SEPT9", "upper"), "SEPTIN9")
  expect_equal(canonicalize_gene_symbols("SEP15", "upper"), "SELENOF")
})

test_that("family casing is normalised to the table convention", {
  # A rodent table carrying MARCHF2 among Marchf1..Marchf11 never pools with
  # the other studies, because the join is a case-sensitive string match.
  expect_equal(canonicalize_gene_symbols("MARCHF2", "title"), "Marchf2")
  expect_equal(canonicalize_gene_symbols("Septin7", "upper"), "SEPTIN7")
})

test_that("symbols outside the two families are untouched", {
  ids <- c("Tnf", "Actb", "10788858", "Rnaset2-ps1", "LOC100038947", "Snord33")
  expect_equal(canonicalize_gene_symbols(ids, "title"), ids)
})

test_that("NA and empty ids survive unchanged", {
  expect_equal(canonicalize_gene_symbols(c(NA, "", "1-Mar"), "title"),
               c(NA, "", "Marchf1"))
})

test_that("a vector is canonicalised elementwise", {
  expect_equal(
    canonicalize_gene_symbols(c("1-Mar", "Sept1", "Tnf", "15-Sep"), "title"),
    c("Marchf1", "Septin1", "Tnf", "Selenof")
  )
})


# ---------------------------------------------------------------------------
# Collision handling
# ---------------------------------------------------------------------------

test_that("canonicalisation can collapse two ids onto one symbol", {
  # A study carrying both 1-Sep and Sept1 ends up with Septin1 twice; the
  # caller must re-aggregate or the study contributes the same gene to a
  # meta-analysis as two independent observations.
  out <- canonicalize_gene_symbols(c("1-Sep", "Sept1"), "title")
  expect_equal(out, c("Septin1", "Septin1"))
  expect_equal(anyDuplicated(out) > 0, TRUE)
})


# ---------------------------------------------------------------------------
# canonicalize_study_symbols — repair plus re-aggregation
# ---------------------------------------------------------------------------

make_study <- function(feature_id, pval) {
  data.frame(feature_id = feature_id, pval = pval,
             effect_size = seq_along(feature_id) / 10,
             stringsAsFactors = FALSE)
}

test_that("a study holding two spellings of one gene keeps a single row", {
  df <- make_study(c("1-Sep", "Sept1", "Tnf"), c(0.4, 0.01, 0.2))
  out <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)

  expect_equal(sum(out$df$feature_id == "Septin1"), 1)
  # The surviving row is the smaller p-value, matching aggregate_to_gene.
  expect_equal(out$df$pval[out$df$feature_id == "Septin1"], 0.01)
  expect_equal(nrow(out$df), 2)
})

test_that("original_feature_id records what the symbol used to be", {
  df <- make_study(c("1-Mar", "Tnf"), c(0.1, 0.2))
  out <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)
  expect_equal(out$df$original_feature_id, c("1-Mar", "Tnf"))
  expect_equal(out$df$feature_id, c("Marchf1", "Tnf"))
})

test_that("the change log names every repaired symbol once", {
  df <- make_study(c("1-Mar", "2-Mar", "1-Mar", "Tnf"), c(0.1, 0.2, 0.3, 0.4))
  out <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)
  expect_equal(nrow(out$changes), 2)
  expect_setequal(out$changes$from, c("1-Mar", "2-Mar"))
  expect_true(all(out$changes$study_id == "GSE1"))
})

test_that("a clean study is returned untouched with no changes", {
  df <- make_study(c("Tnf", "Actb"), c(0.1, 0.2))
  out <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)
  expect_equal(out$df, df)
  expect_equal(nrow(out$changes), 0)
})

test_that("canonicalisation is idempotent", {
  df <- make_study(c("1-Sep", "Sept1", "MARCHF2", "Tnf"), c(0.4, 0.01, 0.2, 0.3))
  once  <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)$df
  twice <- canonicalize_study_symbols(once, "GSE1", verbose = FALSE)$df
  expect_equal(twice$feature_id, once$feature_id)
  expect_equal(nrow(twice), nrow(once))
})

test_that("row order survives the merge", {
  df <- make_study(c("Actb", "1-Sep", "Tnf", "Sept1"), c(0.5, 0.4, 0.2, 0.01))
  out <- canonicalize_study_symbols(df, "GSE1", verbose = FALSE)$df
  expect_equal(out$feature_id, c("Actb", "Tnf", "Septin1"))
})
