library(testthat)

# Source da_utils directly (no installed package needed)
REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "da_utils.R"))


# ---------------------------------------------------------------------------
# assign_groups
# ---------------------------------------------------------------------------

make_pheno <- function(char_values) {
  data.frame(
    characteristics_ch1 = char_values,
    stringsAsFactors = FALSE,
    row.names = paste0("S", seq_along(char_values))
  )
}

test_that("assign_groups labels cases and controls correctly", {
  pheno <- make_pheno(c("chronic pain patient", "healthy control", "pain patient", "normal"))
  cfg_labels <- list(
    case_labels = c("chronic pain", "pain patient"),
    control_labels = c("healthy", "normal")
  )
  out <- assign_groups(pheno, cfg_labels$case_labels, cfg_labels$control_labels)
  expect_equal(out$group[1], "case")
  expect_equal(out$group[2], "control")
  expect_equal(out$group[3], "case")
  expect_equal(out$group[4], "control")
})

test_that("assign_groups returns NA for unmatched samples", {
  pheno <- make_pheno(c("unknown condition", "chronic pain"))
  out <- assign_groups(pheno, "chronic pain", "healthy")
  expect_true(is.na(out$group[1]))
  expect_equal(out$group[2], "case")
})

test_that("assign_groups is case-insensitive", {
  pheno <- make_pheno(c("CHRONIC PAIN PATIENT", "Healthy Control"))
  out <- assign_groups(pheno, "chronic pain", "healthy")
  expect_equal(out$group[1], "case")
  expect_equal(out$group[2], "control")
})


# ---------------------------------------------------------------------------
# build_effect_df
# ---------------------------------------------------------------------------

test_that("build_effect_df produces correct columns", {
  df <- build_effect_df(
    feature_id = c("GENE1", "GENE2"),
    effect_size = c(1.2, -0.5),
    se = c(0.3, 0.1),
    pval = c(0.01, 0.2),
    padj = c(0.05, 0.5),
    n_case = 20L, n_control = 20L,
    study_id = "GSE99999",
    modality = "transcriptomics",
    id_space = "HGNC",
    effect_unit = "log2FC"
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
  expect_true(all(c("feature_id", "effect_size", "se", "pval", "padj",
                    "n_case", "n_control", "study_id", "modality",
                    "id_space", "effect_unit") %in% colnames(df)))
  expect_equal(df$effect_size[1], 1.2)
  expect_equal(df$id_space[1], "HGNC")
})

test_that("build_effect_df coerces types correctly", {
  df <- build_effect_df("GENE1", "1.5", "0.2", "0.05", "0.1",
                        "10", "10", "GSE1", "transcriptomics", "HGNC", "log2FC")
  expect_type(df$effect_size, "double")
  expect_type(df$n_case, "integer")
  expect_type(df$feature_id, "character")
})


# ---------------------------------------------------------------------------
# save_effects / load
# ---------------------------------------------------------------------------

test_that("save_effects writes readable CSV", {
  tmp <- tempdir()
  df <- build_effect_df("GENE1", 1.0, 0.2, 0.01, 0.05,
                        10L, 10L, "GSE1", "transcriptomics", "HGNC", "log2FC")
  save_effects(df, tmp, "transcriptomics", "GSE1")
  out_path <- file.path(tmp, "transcriptomics", "GSE1_effects.csv")
  expect_true(file.exists(out_path))
  loaded <- read.csv(out_path, stringsAsFactors = FALSE)
  expect_equal(loaded$feature_id[1], "GENE1")
  expect_equal(loaded$effect_size[1], 1.0)
})


# ---------------------------------------------------------------------------
# load_included_accessions
#
# Regression cover for plan/2026-08-28-step05-screening-mismatch.md: step 05
# read literature/prisma/screening.csv, whose `uid` column holds PubMed ids,
# and matched it against GEO/PRIDE/MetaboLights accessions. The comparison
# never errored — it just narrowed the corpus to whatever accessions happened
# to have been typed into the literature sheet.
# ---------------------------------------------------------------------------

write_candidates <- function(dir, modality, accessions, verdicts) {
  path <- file.path(dir, paste0(modality, "_candidates.csv"))
  write.csv(
    data.frame(accession = accessions, verdict = verdicts,
               reason_code = "test", stringsAsFactors = FALSE),
    path, row.names = FALSE
  )
  path
}

write_screening_sheet <- function(dir, uids) {
  write.csv(
    data.frame(uid = uids, screening_decision = "include",
               stringsAsFactors = FALSE),
    file.path(dir, "screening.csv"), row.names = FALSE
  )
}

test_that("load_included_accessions returns the modality's included accessions", {
  dir <- withr::local_tempdir()
  write_candidates(dir, "transcriptomics",
                   c("GSE1", "GSE2", "GSE3", "GSE4"),
                   c("include", "exclude", "include", "pending"))
  expect_equal(load_included_accessions(dir, "transcriptomics"),
               c("GSE1", "GSE3"))
})

test_that("load_included_accessions falls back to the manual review file", {
  dir <- withr::local_tempdir()
  write.csv(
    data.frame(accession = c("PXD1", "PXD2"), verdict = c("include", "exclude"),
               stringsAsFactors = FALSE),
    file.path(dir, "proteomics_manual_review.csv"), row.names = FALSE
  )
  expect_equal(load_included_accessions(dir, "proteomics"), "PXD1")
})

test_that("load_included_accessions returns NULL when the modality has no record", {
  dir <- withr::local_tempdir()
  expect_null(load_included_accessions(dir, "genomics"))
})

test_that("load_included_accessions ignores the PubMed literature sheet", {
  # The defect: screening.csv alone must not be treated as a source of
  # accessions. Returning its uids here is what dropped 270 of 283 studies.
  dir <- withr::local_tempdir()
  write_screening_sheet(dir, c("42137963", "41953845", "GSE113941"))
  expect_null(load_included_accessions(dir, "transcriptomics"))
})


# ---------------------------------------------------------------------------
# assert_screening_identifier_space
# ---------------------------------------------------------------------------

test_that("assert_screening_identifier_space accepts a matching space", {
  expect_silent(
    assert_screening_identifier_space(c("GSE1", "GSE3"),
                                      c("GSE1", "GSE2", "GSE3"),
                                      "transcriptomics")
  )
})

test_that("assert_screening_identifier_space rejects a mostly-foreign list", {
  # The transcriptomics shape of the defect: 13 of 47 included ids are GEO
  # accessions, the other 34 are PubMed ids. Prefixes overlap, so a prefix
  # intersection test alone passes; the list is still the wrong file.
  included <- c(paste0("GSE", 1:13), as.character(42137963:42137996))
  expect_error(
    assert_screening_identifier_space(included, paste0("GSE", 1:283),
                                      "transcriptomics"),
    "identifier space"
  )
})

test_that("assert_screening_identifier_space rejects a wholly disjoint space", {
  # The genomics/proteomics/metabolomics shape: nothing can ever match.
  expect_error(
    assert_screening_identifier_space(as.character(42137963:42137996),
                                      paste0("GCST", 1:36), "genomics"),
    "identifier space"
  )
})

test_that("assert_screening_identifier_space rejects an empty intersection", {
  # Right space, but no included accession is actually harmonized: processing
  # zero studies is never the intent when a screening record lists includes.
  expect_error(
    assert_screening_identifier_space(paste0("GSE", 900:910),
                                      paste0("GSE", 1:283), "transcriptomics"),
    "no harmonized record"
  )
})

test_that("assert_screening_identifier_space passes an empty included list through", {
  # NULL means "no screening record for this modality" — not an error.
  expect_silent(
    assert_screening_identifier_space(NULL, paste0("GSE", 1:10), "transcriptomics")
  )
})

test_that("load_included_accessions does not emit NA for a missing verdict", {
  dir <- withr::local_tempdir()
  write.csv(
    data.frame(accession = c("GSE1", "GSE2", "GSE3"),
               verdict = c("include", NA, "include"), stringsAsFactors = FALSE),
    file.path(dir, "transcriptomics_candidates.csv"), row.names = FALSE
  )
  expect_equal(load_included_accessions(dir, "transcriptomics"), c("GSE1", "GSE3"))
})


# ---------------------------------------------------------------------------
# accession_screened_in
# ---------------------------------------------------------------------------

test_that("accession_screened_in keeps an accession named in the record", {
  expect_equal(accession_screened_in(c("GSE1", "GSE2"), c("GSE1")),
               c(TRUE, FALSE))
})

test_that("accession_screened_in admits a derived unit via its parent", {
  # Only parents are ever screened; a regional split appears in no sheet.
  included <- c("GSE180627", "GSE197233")
  units <- c("GSE180627_PAG", "GSE180627_S1", "GSE197233_mPFC")
  expect_true(all(accession_screened_in(units, included)))
})

test_that("accession_screened_in does not admit a derived unit of an excluded parent", {
  expect_false(accession_screened_in("GSE999_PAG", c("GSE180627")))
})

test_that("accession_screened_in resolves a twice-derived unit", {
  expect_true(accession_screened_in("GSE1_DRG_left", c("GSE1")))
})

test_that("accession_screened_in returns a logical of the input length", {
  expect_length(accession_screened_in(character(0), c("GSE1")), 0)
  expect_type(accession_screened_in("GSE1", c("GSE1")), "logical")
})

# ---------------------------------------------------------------------------
# assign_groups: a control arm that names the disease it does not have
# ---------------------------------------------------------------------------
# Matching is fixed substring and first-match-wins, so once "endometriosis" is
# a case label the value "no endometriosis" hits it first and the control
# becomes a case. The study then has two case arms and no split, which is what
# hid a body of endometriosis studies behind
# `no_case_control_split_in_series_matrix`.

CASE_LBL <- c("chronic pain", "endometriosis")
CTRL_LBL <- c("healthy", "control", "normal")

pheno_from <- function(values) {
  data.frame(characteristics_ch1 = values, stringsAsFactors = FALSE)
}

test_that("a negated disease label is a control, not a case", {
  out <- assign_groups(
    pheno_from(c("disease: endometriosis", "disease: no endometriosis",
                 "disease: non-endometriosis", "disease: healthy")),
    CASE_LBL, CTRL_LBL)
  expect_equal(out$group, c("case", "control", "control", "control"))
})

test_that("every negation form the deposits use is recognised", {
  out <- assign_groups(
    pheno_from(c("non endometriosis", "non_endometriosis",
                 "without endometriosis", "endometriosis-free",
                 "disease free endometrium", "no pain")),
    CASE_LBL, CTRL_LBL)
  expect_true(all(out$group == "control"))
})

test_that("the disease itself is still a case", {
  out <- assign_groups(pheno_from(c("endometriosis", "chronic pain")),
                       CASE_LBL, CTRL_LBL)
  expect_true(all(out$group == "case"))
})

test_that("a pre-assigned group column is still returned untouched", {
  # Every included study carries a curated groups.csv, and that path must not
  # be reachable by the vocabulary at all.
  pheno <- data.frame(group = c("case", "control"), stringsAsFactors = FALSE)
  expect_equal(assign_groups(pheno, CASE_LBL, CTRL_LBL)$group,
               c("case", "control"))
})
