library(testthat)

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
source(file.path(REPO_ROOT, "R", "meta_utils.R"))


# ---------------------------------------------------------------------------
# adjust_pooled_fdr
# ---------------------------------------------------------------------------

test_that("adjust_pooled_fdr matches p.adjust on the k >= min_k subset alone", {
  df <- data.frame(
    feature_id = paste0("G", 1:6),
    k          = c(3L, 2L, 1L, 4L, 1L, 2L),
    pval       = c(1e-6, 1e-3, 1e-9, 0.02, 1e-8, 0.5)
  )
  out <- adjust_pooled_fdr(df)
  keep <- df$k >= 2
  expect_equal(out$padj[keep], p.adjust(df$pval[keep], method = "BH"))
})

test_that("adjust_pooled_fdr leaves single-study rows at NA", {
  df <- data.frame(k = c(1L, 2L, 1L), pval = c(1e-9, 0.04, 1e-8))
  out <- adjust_pooled_fdr(df)
  expect_true(all(is.na(out$padj[df$k == 1])))
  expect_false(is.na(out$padj[df$k == 2]))
})

test_that("single-study rows do not change the adjusted p-values of pooled rows", {
  pooled_only <- data.frame(k = c(2L, 3L, 2L), pval = c(1e-5, 1e-3, 0.04))
  with_k1 <- rbind(
    pooled_only,
    data.frame(k = rep(1L, 500), pval = 10^-runif(500, 1, 12))
  )
  expect_equal(
    adjust_pooled_fdr(with_k1)$padj[1:3],
    adjust_pooled_fdr(pooled_only)$padj
  )
})

test_that("adjust_pooled_fdr preserves row order and adds no rows", {
  df <- data.frame(
    feature_id = c("Zzz", "Aaa", "Mmm"),
    k          = c(2L, 1L, 3L),
    pval       = c(0.3, 0.001, 0.002)
  )
  out <- adjust_pooled_fdr(df)
  expect_equal(nrow(out), 3)
  expect_equal(out$feature_id, c("Zzz", "Aaa", "Mmm"))
})

test_that("adjust_pooled_fdr excludes NA p-values from the family", {
  df <- data.frame(k = c(2L, 2L, 2L), pval = c(0.01, NA, 0.02))
  out <- adjust_pooled_fdr(df)
  expect_true(is.na(out$padj[2]))
  expect_equal(out$padj[c(1, 3)],
               p.adjust(c(0.01, 0.02), method = "BH"))
})

test_that("adjust_pooled_fdr returns all-NA when no row qualifies", {
  df <- data.frame(k = c(1L, 1L), pval = c(1e-9, 0.5))
  out <- adjust_pooled_fdr(df)
  expect_true(all(is.na(out$padj)))
})

test_that("adjust_pooled_fdr handles a zero-row frame without error", {
  df <- data.frame(k = integer(0), pval = numeric(0))
  out <- adjust_pooled_fdr(df)
  expect_equal(nrow(out), 0)
  expect_true("padj" %in% colnames(out))
})

test_that("adjust_pooled_fdr honours custom column names and min_k", {
  df <- data.frame(kk = c(3L, 2L, 3L), p = c(0.01, 0.2, 0.3))
  out <- adjust_pooled_fdr(df, pval_col = "p", k_col = "kk", min_k = 3L,
                           out_col = "q")
  expect_true(is.na(out$q[2]))
  expect_equal(out$q[c(1, 3)], p.adjust(c(0.01, 0.3), method = "BH"))
})


# ---------------------------------------------------------------------------
# usable_se / Z_975
# ---------------------------------------------------------------------------

test_that("usable_se keeps only rows with a positive, non-missing SE", {
  df <- data.frame(id = 1:5, se = c(0.2, NA, 0, -0.1, 1e-8))
  expect_equal(usable_se(df)$id, c(1L, 5L))
})

test_that("usable_se honours a custom column name", {
  df <- data.frame(id = 1:3, sigma = c(0.1, NA, -2))
  expect_equal(usable_se(df, se_col = "sigma")$id, 1L)
})

test_that("usable_se returns a zero-row frame rather than erroring", {
  df <- data.frame(id = integer(0), se = numeric(0))
  expect_equal(nrow(usable_se(df)), 0)
})

test_that("Z_975 is qnorm(0.975) rounded exactly as the Python side rounds it", {
  # harmonize.py hardcodes 1.959964; both sides must agree to that precision,
  # which is finer than any effect-size interval the pipeline reports.
  expect_equal(Z_975, qnorm(0.975), tolerance = 1e-6)
  expect_identical(Z_975, 1.959964)
})


# ---------------------------------------------------------------------------
# moderator_df
# ---------------------------------------------------------------------------

test_that("moderator_df returns the moderator df, not the coefficient count", {
  skip_if_not_installed("metafor")
  suppressPackageStartupMessages(library(metafor))
  set.seed(1)
  y <- rnorm(12); s <- runif(12, .1, .3)
  g <- factor(rep(c("a", "b", "c"), 4))
  res <- rma(yi = y, sei = s, mods = ~ g, method = "REML")
  expect_equal(moderator_df(res), 2)          # 3 groups -> 2 df
  expect_equal(moderator_df(res), unname(res$QMdf[1]))
  expect_false(moderator_df(res) == res$p)    # res$p is 3 (intercept + 2)
})

test_that("moderator_df is unaffected by unused factor levels", {
  skip_if_not_installed("metafor")
  suppressPackageStartupMessages(library(metafor))
  set.seed(2)
  y <- rnorm(8); s <- runif(8, .1, .3)
  g <- droplevels(factor(rep(c("a", "b"), 4), levels = c("a", "b", "c", "d")))
  res <- rma(yi = y, sei = s, mods = ~ g, method = "REML")
  expect_equal(moderator_df(res), 1)
})


# ---------------------------------------------------------------------------
# top_pooled_features
# ---------------------------------------------------------------------------

test_that("top_pooled_features ranks within the adjusted family only", {
  df <- data.frame(
    feature_id  = c("k1a", "p1", "k1b", "p2", "p3"),
    pval_pooled = c(1e-30, 1e-10, 1e-20, 1e-8, 1e-6),
    padj_pooled = c(NA, 1e-6, NA, 1e-5, 1e-3)   # k=1 rows carry NA
  )
  # k1a/k1b have the smallest p-values but never entered a meta-analysis
  expect_equal(top_pooled_features(df, 3), c("p1", "p2", "p3"))
})

test_that("top_pooled_features returns at most n, ordered by p-value", {
  df <- data.frame(
    feature_id  = c("b", "a", "c"),
    pval_pooled = c(0.02, 0.001, 0.3),
    padj_pooled = c(0.05, 0.01, 0.4)
  )
  expect_equal(top_pooled_features(df, 2), c("a", "b"))
  expect_equal(top_pooled_features(df, 99), c("a", "b", "c"))
})

test_that("top_pooled_features returns an empty vector when nothing pooled", {
  df <- data.frame(feature_id = c("x", "y"),
                   pval_pooled = c(1e-9, 1e-8),
                   padj_pooled = c(NA_real_, NA_real_))
  expect_equal(top_pooled_features(df, 5), character(0))
})

test_that("top_pooled_features handles a zero-row frame", {
  df <- data.frame(feature_id = character(0), pval_pooled = numeric(0),
                   padj_pooled = numeric(0))
  expect_equal(top_pooled_features(df, 5), character(0))
})


# ---------------------------------------------------------------------------
# superseded_studies
# ---------------------------------------------------------------------------

test_that("superseded_studies returns the excluded accessions", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    study_id = c("GSE1", "GSE2"),
    superseded_by = c("GSE1_A;GSE1_B", "GSE2_A"),
    reason = c("split by tissue", "split by region")
  ), tmp, row.names = FALSE)
  expect_equal(superseded_studies(tmp), c("GSE1", "GSE2"))
})

test_that("superseded_studies returns empty when the registry is absent", {
  expect_equal(superseded_studies(tempfile(fileext = ".csv")), character(0))
})

test_that("superseded_studies returns empty for a header-only registry", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(study_id = character(0), superseded_by = character(0),
                       reason = character(0)), tmp, row.names = FALSE)
  expect_equal(superseded_studies(tmp), character(0))
})

# --- study_strata.csv, the config-driven stratum map -------------------------
# Stratum assignment used to be a hardcoded vector in 06b. A search amendment
# added 46 studies at once, which would have meant editing analysis source to
# record a retrieval decision.

test_that("study_strata.csv is well formed and every stratum is named", {
  path <- file.path(REPO_ROOT, "conf", "analysis", "study_strata.csv")
  skip_if_not(file.exists(path), "study_strata.csv absent")
  sf <- read.csv(path, stringsAsFactors = FALSE)
  expect_true(all(c("study_id", "stratum") %in% names(sf)))
  expect_true(all(nzchar(sf$study_id)))
  expect_true(all(nzchar(sf$stratum)))
  # One row per study: a study in two strata would enter the pool twice.
  expect_equal(anyDuplicated(sf$study_id), 0)
})

# Coverage, not just shape. A study admitted by screening and absent from both
# stratum maps is dropped from every model stratum by 06b, which reports it and
# carries on -- so the stratum count stays where it was while the corpus grows,
# and nothing fails. That is how GSE276193 was admitted to the transcriptomic
# arm on 2026-09-09 and still contributed to no stratum: it belongs to
# endometriosis_human alongside GSE205494 and 25 others and was in neither map.
test_that("every included transcriptomic study is assigned a stratum", {
  strata <- file.path(REPO_ROOT, "conf", "analysis", "study_strata.csv")
  cands  <- file.path(REPO_ROOT, "literature", "prisma",
                      "transcriptomics_candidates.csv")
  sup_p  <- file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv")
  uns_p  <- file.path(REPO_ROOT, "conf", "analysis", "unstratified_studies.csv")
  skip_if_not(all(file.exists(c(strata, cands, sup_p, uns_p))))

  # STUDY_MODEL is a vector inside 06b, which cannot be sourced here (it loads
  # metafor and dplyr at the top level), so lift the literal out of the source.
  src <- readLines(file.path(REPO_ROOT, "pipeline", "06b_stratified_meta.R"),
                   warn = FALSE)
  at <- grep("^STUDY_MODEL <- c\\(", src)
  expect_length(at, 1)
  close_at <- grep("^\\)", src); close_at <- close_at[close_at > at][1]
  env <- new.env()
  eval(parse(text = paste(src[at:close_at], collapse = "\n")), envir = env)

  mapped <- union(names(env$STUDY_MODEL),
                  read.csv(strata, stringsAsFactors = FALSE)$study_id)
  included <- {
    cd <- read.csv(cands, stringsAsFactors = FALSE)
    cd$accession[cd$verdict %in% "include"]
  }
  superseded <- read.csv(sup_p, stringsAsFactors = FALSE)$study_id

  # Studies deliberately left out of every stratum, each with a reason. The
  # registry exists so that "in no stratum" is a decision on the record rather
  # than an omission, which is the same distinction superseded_studies.csv
  # draws. It is subtracted here and nowhere else: 06b still reports these
  # studies, because a reader of its output should see the gap.
  deliberate <- read.csv(uns_p, stringsAsFactors = FALSE)$study_id

  unstratified <- sort(setdiff(setdiff(setdiff(included, superseded), mapped),
                               deliberate))
  expect_equal(unstratified, character(0),
               info = paste("admitted by screening but in no stratum map, so",
                            "06b will drop them from every stratum:",
                            paste(unstratified, collapse = ", ")))

  # The registry must not outlive its entries. A study that leaves the corpus,
  # gains a stratum, or is superseded should drop out of it, or the exemption
  # silently covers a case nobody has looked at since.
  stale <- sort(intersect(deliberate, union(mapped, superseded)))
  expect_equal(stale, character(0),
               info = paste("listed in unstratified_studies.csv but now",
                            "stratified or superseded; remove them:",
                            paste(stale, collapse = ", ")))
  absent <- sort(setdiff(deliberate, included))
  expect_equal(absent, character(0),
               info = paste("listed in unstratified_studies.csv but not",
                            "included in the corpus:",
                            paste(absent, collapse = ", ")))
})

test_that("no study is both superseded and assigned a stratum", {
  strata <- file.path(REPO_ROOT, "conf", "analysis", "study_strata.csv")
  sup <- file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv")
  skip_if_not(file.exists(strata) && file.exists(sup))
  sf <- read.csv(strata, stringsAsFactors = FALSE)
  ss <- read.csv(sup, stringsAsFactors = FALSE)
  # A superseded study given a stratum would be excluded by the loader while
  # appearing in the stratified design, which reads as a silent inconsistency.
  expect_length(intersect(sf$study_id, ss$study_id), 0)
})


# ---------------------------------------------------------------------------
# conf/analysis/species_overrides.csv
#
# 06g assigns species from the manifest, which records what GEO says about a
# series. GSE236754 is listed as Mus musculus + Rattus norvegicus over 56
# samples, so it was excluded as multi-species -- but only 32 samples carry a
# case/control assignment, 05d filters to those, and the matrix analysed is
# entirely rat. It was withheld from the pool it belongs to. The override
# records what was analysed without editing the manifest.
# ---------------------------------------------------------------------------

test_that("species_overrides.csv has the columns 06g reads", {
  p <- file.path(REPO_ROOT, "conf", "analysis", "species_overrides.csv")
  skip_if_not(file.exists(p))
  ovr <- read.csv(p, stringsAsFactors = FALSE)
  expect_true(all(c("study_id", "species", "reason") %in% colnames(ovr)))
  expect_true(all(nzchar(ovr$study_id)), info = "no blank study_id")
  expect_true(all(nzchar(ovr$species)), info = "no blank species")
  # A reason is the point: an override with no justification is indistinguishable
  # from a typo, and this file silently changes which pool a study joins.
  expect_true(all(nchar(ovr$reason) > 40), info = "every override needs a reason")
})

test_that("06g applies the overrides before counting species", {
  src <- readLines(file.path(REPO_ROOT, "pipeline", "06g_species_meta.R"), warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]
  ovr <- grep("species_overrides.csv", code)
  nsp <- grep("n_species\\s*<- vapply", code)
  expect_gt(length(ovr), 0)
  expect_length(nsp, 1)
  expect_lt(min(ovr), nsp)
})


# ---------------------------------------------------------------------------
# drop_sparse_for_clustering
#
# 06b's cross-model heatmap calls heatmap(), which clusters rows and columns
# with dist(). A pair of rows sharing no non-NA column yields NaN, and hclust
# then aborts with "NA/NaN/Inf in foreign function call (arg 10)" -- which is
# what killed the 2026-08-29 run after every stratum had already been pooled.
# The matrix went sparse because the search amendment took the design from 9
# strata to 14, several of them human-only, so a top feature is absent from
# most columns.
# ---------------------------------------------------------------------------

test_that("drop_sparse_for_clustering keeps a dense matrix intact", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3,
              dimnames = list(c("a", "b", "c"), c("x", "y")))
  expect_equal(drop_sparse_for_clustering(m), m)
})

test_that("drop_sparse_for_clustering removes an all-NA column", {
  m <- matrix(c(1, 2, NA, NA, 3, 4), nrow = 2,
              dimnames = list(c("a", "b"), c("x", "y", "z")))
  got <- drop_sparse_for_clustering(m)
  expect_equal(colnames(got), c("x", "z"))
})

test_that("drop_sparse_for_clustering removes a row below the observation floor", {
  m <- matrix(c(1, NA, 2, 3, 4, 5), nrow = 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("x", "y")))
  got <- drop_sparse_for_clustering(m, min_obs = 2)
  expect_equal(rownames(got), c("b", "c"))
})

test_that("drop_sparse_for_clustering declines when one row survives the floor", {
  # Trimming can leave a matrix too small to cluster at all; that is a decline,
  # not a one-row heatmap.
  m <- matrix(c(1, NA, 2, 3, NA, 4), nrow = 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("x", "y")))
  expect_null(drop_sparse_for_clustering(m, min_obs = 2))
})

test_that("drop_sparse_for_clustering returns NULL rather than an unclusterable matrix", {
  # Two rows that share no observed column: dist() gives NaN whatever is
  # dropped, so the caller must be told to skip clustering, not handed a
  # matrix that will abort inside hclust.
  m <- matrix(c(1, NA, NA, 2), nrow = 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("x", "y")))
  expect_null(drop_sparse_for_clustering(m, min_obs = 1))
})

test_that("drop_sparse_for_clustering returns NULL when too little survives", {
  m <- matrix(c(1, NA, NA, NA), nrow = 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("x", "y")))
  expect_null(drop_sparse_for_clustering(m))
})


# ---------------------------------------------------------------------------
# resolve_by_parent
# ---------------------------------------------------------------------------

test_that("resolve_by_parent inherits a derived unit's value from its parent", {
  map <- c(GSE180627 = "CCI", GSE197233 = "SNI")
  expect_equal(
    resolve_by_parent(c("GSE180627_PAG", "GSE197233_mPFC", "GSE180627"), map),
    c("CCI", "SNI", "CCI")
  )
})

test_that("resolve_by_parent accepts a region token carrying a digit", {
  # GSE180627_S1 is the primary somatosensory cortex. A rule requiring the
  # token to be purely alphabetic left this one unit of four unresolved.
  map <- c(GSE180627 = "CCI")
  expect_equal(resolve_by_parent("GSE180627_S1", map), "CCI")
})

test_that("resolve_by_parent prefers an exact hit to an inherited one", {
  # A split that changes the mapped property is recorded per unit, and must
  # not be overridden by the parent.
  map <- c(GSE241361 = "SNI", GSE241361_DRG = "CCI")
  expect_equal(resolve_by_parent("GSE241361_DRG", map), "CCI")
})

test_that("resolve_by_parent strips repeated suffixes", {
  map <- c(GSE1 = "SNI")
  expect_equal(resolve_by_parent("GSE1_DRG_left", map), "SNI")
})

test_that("resolve_by_parent returns NA for an unmapped study", {
  map <- c(GSE180627 = "CCI")
  expect_true(is.na(resolve_by_parent("GSE999999_PAG", map)))
  expect_true(is.na(resolve_by_parent("GSE999999", map)))
})

test_that("resolve_by_parent handles an empty input without looping", {
  expect_length(resolve_by_parent(character(0), c(GSE1 = "SNI")), 0)
})


# The forest and the funnel of a modality are a pair, and used to describe
# different genes: the forest drew top_pooled_features()[1] while the funnel
# took the most-measured feature, so the single-cell arm's forest showed Ppm1l
# and its funnel showed C3. Nothing in either filename or title said so.
test_that("the funnel and the forest select the same feature", {
  src <- readLines(file.path(REPO_ROOT, "pipeline", "06_meta_analysis.R"),
                   warn = FALSE)
  at <- grep("^plot_funnel <- function", src)
  expect_length(at, 1)
  close_at <- grep("^\\}", src); close_at <- close_at[close_at > at][1]
  body <- src[at:close_at]
  expect_true(any(grepl("top_pooled_features\\(pooled_df", body)))
  # The old selection ranked by how many studies measured a feature, which is
  # not a ranking within the adjusted family and picks a different gene.
  expect_false(any(grepl("count\\(feature_id", body)))
  # Both take the pooled table, so neither can rank across unadjusted rows.
  expect_true(any(grepl("plot_funnel\\(effects_df, pooled,", src)))
})
