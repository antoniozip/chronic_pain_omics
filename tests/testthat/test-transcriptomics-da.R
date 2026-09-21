library(testthat)

# transcriptomics_da.R cannot be source()d here: its top-level block loads
# GEOquery/DESeq2/edgeR, which are absent under R 4.3.3 (renv.lock pins
# Bioconductor 3.22, requiring R >= 4.5). These are source-level regression
# guards over the SE conventions, which are the only thing a wrong edit
# could silently break. See plan/cohort_expansion_plan.md §7.5.

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
SRC <- readLines(file.path(REPO_ROOT, "R", "transcriptomics_da.R"), warn = FALSE)
CODE <- SRC[!grepl("^\\s*#", SRC)]   # strip comment lines


test_that("apeglm shrinkage never returns", {
  # lfcShrink(type='apeglm') yields a posterior SD of ~0.001 against a true SE
  # of 0.10-0.25, which would collapse random-effects weights onto a few studies.
  expect_false(any(grepl("lfcShrink", CODE)))
  expect_false(any(grepl("apeglm", CODE)))
})

test_that("DESeq2 SE is taken from the unshrunken results() object", {
  expect_true(any(grepl("results\\(dds", CODE)))
})

test_that("limma SE delegates to se_from_t, not an inline qt formula", {
  expect_false(any(grepl("qt\\(0\\.975", CODE)))
  expect_true(any(grepl("se_from_t\\(", CODE)))
})

test_that("the log2 heuristic delegates to maybe_log2", {
  expect_true(any(grepl("maybe_log2\\(", CODE)))
  # The bare threshold must not be re-hardcoded at 100.
  expect_false(any(grepl("q75 > 100", CODE)))
})

test_that("se_utils and cuffdiff_da are sourced", {
  expect_true(any(grepl("se_utils\\.R", SRC)))
  expect_true(any(grepl("cuffdiff_da\\.R", SRC)))
})

test_that("the header docstring does not advertise apeglm shrinkage", {
  header <- paste(SRC[1:10], collapse = " ")
  # The header may NAME apeglm in order to warn against it; what it must never
  # do is describe apeglm shrinkage as the method in use.
  expect_false(grepl("Apeglm LFC shrinkage", header))
  expect_true(grepl("UNSHRUNKEN", header))
})

test_that("no inline limma SE formula remains", {
  # Presence of se_from_t() is not enough: a partial migration that left one of
  # the four call sites as abs(tt$logFC / tt$t) would otherwise pass silently.
  expect_false(any(grepl("abs\\(tt\\$logFC\\s*/\\s*tt\\$t\\)", CODE)))
})

test_that("the inline q75 log2 threshold block is fully removed", {
  expect_false(any(grepl("\\bq75\\b", CODE)))
})


# ---------------------------------------------------------------------------
# The series-matrix group gate must not pre-empt supplementary count files.
#
# Ten of the 57 selected studies (GSE102721, GSE102937, GSE117526, GSE138024,
# GSE160543, GSE175760, GSE224814, GSE236754, GSE245768, GSE272517) supply
# expression as {GSE}_raw_count_merged.txt.gz, whose columns are named
# naive_control_1 / injured_case_1 rather than GSM ids. run_deseq2_from_matrix()
# groups those from the column names via classify_cols_by_vocab(), consulting
# the series-matrix phenotype only as a fallback. Gating on that phenotype
# first dropped all ten with "insufficient group sizes (case=0, ctrl=0)" and,
# because save_effects() then never fires, left their stale tables in place.
# ---------------------------------------------------------------------------

gate_line <- grep("Skipping: insufficient group sizes", CODE)
suppl_line <- grep("suppl_file <- find_suppl_count_file", CODE)

test_that("the expression source is determined before the group-size gate", {
  expect_length(gate_line, 1)
  expect_length(suppl_line, 1)
  expect_lt(suppl_line, gate_line)
})

test_that("the group-size gate is skipped when a supplementary matrix exists", {
  guard <- grep("is.null(suppl_file)", CODE, fixed = TRUE)
  expect_true(any(guard < gate_line & guard > suppl_line),
              info = "the gate must be guarded by is.null(suppl_file)")
})

test_that("classify_cols_by_vocab groups the count-file column names", {
  # The premise of the fix: once the gate lets these studies through, the
  # downstream path assigns their groups correctly with no annotation file.
  # Evaluate the vocabulary block and the function together: the whole file
  # cannot be sourced here (GEOquery/DESeq2 are absent under R 4.3.3).
  start <- grep("^PAIN_TOKENS_R <- c\\(", SRC)
  fn_at <- grep("^classify_cols_by_vocab <- function", SRC)
  expect_length(start, 1)
  expect_length(fn_at, 1)
  ends <- grep("^\\}", SRC)
  block <- SRC[start:min(ends[ends > fn_at])]
  eval(parse(text = paste(block, collapse = "\n")), envir = environment())
  got <- classify_cols_by_vocab(c("naive_control_1", "naive_control_2",
                                  "injured_case_1", "injured_case_2",
                                  "vehicle_control_1", "cfa_case_1",
                                  "sham_control_1", "cci_case_1"))
  expect_equal(unname(got),
               c("control", "control", "case", "case",
                 "control", "case", "control", "case"))
})


# ---------------------------------------------------------------------------
# Matrix columns keyed by GSM must be joined to the phenotype directly.
#
# The mirror image of the gate bug above. 05d writes {GSE}_raw_count_merged.txt.gz
# with GSM column names for some studies and descriptive names for others.
# run_deseq2_from_matrix() had two strategies -- classify_cols_by_vocab() over
# the column names, and a substring match against pheno titles -- and neither
# can read a GSM. GSE154816, GSE241361, GSE256472 and GSE295863 therefore
# reported "Cannot assign groups (vocab: 0/0; title match: 0/0)" while holding a
# perfectly good GSM-keyed annotation (e.g. GSE256472: 3 case, 3 control over
# 6 GSM columns). The join was available and simply never attempted.
# ---------------------------------------------------------------------------

# The whole file cannot be source()d here (GEOquery/DESeq2 absent under R
# 4.3.3), so lift one function plus the module constants it closes over.
eval_helper <- function(name, consts = character()) {
  ends <- grep("^\\}", SRC)
  take <- function(at) {
    stopifnot(length(at) == 1)
    close_at <- ends[ends > at]
    lines <- if (length(close_at) == 0) at else at:min(close_at)
    SRC[lines]
  }
  text <- character()
  for (cn in consts) {
    text <- c(text, take(grep(paste0("^", cn, " *<- "), SRC)))
  }
  text <- c(text, take(grep(paste0("^", name, " <- function"), SRC)))
  eval(parse(text = paste(text, collapse = "\n")), envir = parent.frame())
}

test_that("groups_from_sample_ids joins GSM columns to phenotype rownames", {
  eval_helper("groups_from_sample_ids")
  pheno <- data.frame(group = c("control", "control", "case", "case"),
                      row.names = c("GSM1", "GSM2", "GSM3", "GSM4"),
                      stringsAsFactors = FALSE)
  got <- groups_from_sample_ids(c("GSM1", "GSM3", "GSM4", "GSM2"), pheno)
  expect_equal(unname(got), c("control", "case", "case", "control"))
})

test_that("groups_from_sample_ids returns NA for columns absent from the phenotype", {
  eval_helper("groups_from_sample_ids")
  pheno <- data.frame(group = c("control", "case"),
                      row.names = c("GSM1", "GSM2"), stringsAsFactors = FALSE)
  got <- groups_from_sample_ids(c("GSM1", "GSM9"), pheno)
  expect_equal(unname(got), c("control", NA_character_))
})

test_that("groups_from_sample_ids yields nothing for descriptive column names", {
  # The ten count-file studies: their columns are naive_control_1 etc. and the
  # phenotype is GSM-keyed, so this strategy must decline and let the vocabulary
  # path handle them rather than mislabelling anything.
  eval_helper("groups_from_sample_ids")
  pheno <- data.frame(group = c("control", "case"),
                      row.names = c("GSM1", "GSM2"), stringsAsFactors = FALSE)
  got <- groups_from_sample_ids(c("naive_control_1", "injured_case_1"), pheno)
  expect_true(all(is.na(got)))
})

test_that("run_deseq2_from_matrix attempts the sample-id join", {
  body_at <- grep("^run_deseq2_from_matrix <- function", CODE)
  vocab_at <- grep("classify_cols_by_vocab(colnames(mat_num))", CODE, fixed = TRUE)
  join_at <- grep("groups_from_sample_ids(colnames(mat_num)", CODE, fixed = TRUE)
  expect_length(join_at, 1)
  expect_gt(join_at, body_at)
  expect_lt(join_at, vocab_at)
})


# ---------------------------------------------------------------------------
# find_suppl_count_file selection
#
# GSE135080 ships both GSE135080_matrix_counts.txt.gz -- whitespace-aligned,
# with a blank leading header cell, and the only file carrying gene symbols --
# and 05d's GSE135080_raw_count_merged.txt.gz, which parses cleanly but has no
# gene-id column at all. Selection stays name-neutral; the whitespace fallback
# above is what makes the symbol-bearing file readable.
# ---------------------------------------------------------------------------

test_that("find_suppl_count_file still returns a lone non-merged candidate", {
  eval_helper("find_suppl_count_file", c("COUNT_PATTERNS_R", "SKIP_PATTERNS_R"))
  dir <- withr::local_tempdir()
  file.create(file.path(dir, "GSE1_matrix_counts.txt.gz"))
  expect_equal(basename(find_suppl_count_file(dir)), "GSE1_matrix_counts.txt.gz")
})


# ---------------------------------------------------------------------------
# Whitespace-aligned count matrices
#
# GSE135080 ships GSE135080_matrix_counts.txt.gz, which is space-aligned with a
# blank leading header cell -- the only file for that study carrying gene
# symbols (Adcy5, Prox2). read.delim() splits on tabs, so it read as a single
# column and the study failed with "Row/ID length mismatch: 0 vs 3858". 05d's
# sibling _raw_count_merged.txt.gz parses cleanly but has no gene-id column at
# all, so its first *count* column became the feature ids (39, 3, 24 ...),
# which is what reached the meta-analysis.
# ---------------------------------------------------------------------------

test_that("the matrix reader falls back to whitespace separation", {
  fallback <- grep("sep = \"\"", CODE, fixed = TRUE)
  expect_gt(length(fallback), 0)
})

test_that("find_suppl_count_file does not prefer the merged matrix by name", {
  # Reverted: 05d's merge can drop the gene-id column entirely, so preferring
  # it by name picks the file *without* symbols. Selection stays name-neutral.
  expect_length(grep("_raw_count_merged", CODE), 0)
})


# ---------------------------------------------------------------------------
# Supplementary-matrix columns joined through per-sample metadata.
#
# A submitter's count matrix names its columns after the library or tube --
# EDCV009_Th1-17_S5, 7780.EE.0002, G1_42L, HC-179_FPKM -- which matches neither
# the GSM accession nor the pain vocabulary. GEO carries that same string in
# !Sample_description, sometimes labelled ("Library name: G1_42L", "Column name
# in counts.txt: 7780.EE.0002"), so the column can be joined to a sample and
# the group taken from the curated annotation instead of being guessed from
# what the label appears to mean. 14 of the 46 studies the 2026-08-28 search
# amendment contributed failed only on this; see
# plan/2026-08-29-amendment-da-run.md.
# ---------------------------------------------------------------------------

meta_pheno <- function(desc, groups = c("control", "control", "case", "case")) {
  data.frame(description = desc, group = groups,
             row.names = paste0("GSM", seq_along(desc)),
             stringsAsFactors = FALSE)
}

test_that("groups_from_sample_metadata matches a description verbatim", {
  eval_helper("groups_from_sample_metadata")
  pheno <- meta_pheno(c("1", "11", "12", "4"))
  got <- groups_from_sample_metadata(c("12", "1", "4", "11"), pheno)
  expect_equal(unname(got), c("case", "control", "case", "control"))
})

test_that("groups_from_sample_metadata strips a labelled prefix", {
  eval_helper("groups_from_sample_metadata")
  pheno <- meta_pheno(c("Library name: G1_42L", "Library name: G1_46L",
                        "Library name: G2_45R", "Library name: G2_49R"))
  got <- groups_from_sample_metadata(c("G1_42L", "G2_45R", "G2_49R", "G1_46L"), pheno)
  expect_equal(unname(got), c("control", "case", "case", "control"))
})

test_that("groups_from_sample_metadata matches a column that extends the label", {
  # GSE333494: description "HC-179", column "HC-179_FPKM".
  eval_helper("groups_from_sample_metadata")
  pheno <- meta_pheno(c("Library name: HC-179", "Library name: HC-183",
                        "Library name: HC-186", "Library name: HC-190"))
  got <- groups_from_sample_metadata(
    c("HC-179_FPKM", "HC-183_FPKM", "HC-186_FPKM", "HC-190_FPKM"), pheno)
  expect_equal(unname(got), c("control", "control", "case", "case"))
})

test_that("groups_from_sample_metadata declines an ambiguous containment", {
  # "E01" is a prefix of two samples: assigning either would be a coin toss.
  eval_helper("groups_from_sample_metadata")
  pheno <- meta_pheno(c("E01_0", "E01_2", "N_D0", "N_D2"))
  got <- groups_from_sample_metadata(c("E01", "N_D0"), pheno)
  expect_equal(unname(got), c(NA_character_, "case"))
})

test_that("groups_from_sample_metadata prefers an exact match over a containment", {
  eval_helper("groups_from_sample_metadata")
  pheno <- meta_pheno(c("S1", "S1_rep2", "S2", "S2_rep2"))
  got <- groups_from_sample_metadata(c("S1", "S2"), pheno)
  expect_equal(unname(got), c("control", "case"))
})

test_that("groups_from_sample_metadata yields nothing without a usable field", {
  eval_helper("groups_from_sample_metadata")
  pheno <- data.frame(group = c("control", "case"), row.names = c("GSM1", "GSM2"),
                      stringsAsFactors = FALSE)
  got <- groups_from_sample_metadata(c("colA", "colB"), pheno)
  expect_true(all(is.na(got)))
})

test_that("groups_from_sample_metadata searches every metadata column present", {
  eval_helper("groups_from_sample_metadata")
  pheno <- data.frame(
    title = c("Normal replicate 1", "Normal replicate 2",
              "Endometriosis replicate 1", "Endometriosis replicate 2"),
    supplementary_file = c("ftp://x/GSM1_ctrlA.txt.gz", "ftp://x/GSM2_ctrlB.txt.gz",
                           "ftp://x/GSM3_caseA.txt.gz", "ftp://x/GSM4_caseB.txt.gz"),
    group = c("control", "control", "case", "case"),
    row.names = paste0("GSM", 1:4), stringsAsFactors = FALSE)
  got <- groups_from_sample_metadata(c("ctrlA", "caseB"), pheno)
  expect_equal(unname(got), c("control", "case"))
})

test_that("run_deseq2_from_matrix tries the metadata join after the title match", {
  body_at  <- grep("^run_deseq2_from_matrix <- function", CODE)
  title_at <- grep("Title match", CODE, fixed = TRUE)
  meta_at  <- grep("groups_from_sample_metadata(colnames(mat_num)", CODE, fixed = TRUE)
  expect_length(meta_at, 1)
  expect_gt(meta_at, body_at)
  expect_gt(meta_at, min(title_at))
})


test_that("the numeric conversion keeps the matrix column names verbatim", {
  # as.data.frame() runs make.names() unless told not to, turning GSE343056's
  # "con-1" into "con.1" and GSE333494's "HC-186_Raw.Read.Count" into
  # "HC.186_Raw.Read.Count". Every column-to-sample join then compares a
  # mangled name against GEO's unmangled one and finds nothing, which is
  # indistinguishable from a study whose metadata simply does not match.
  at <- grep("mat_num <- suppressWarnings(", CODE, fixed = TRUE)
  expect_length(at, 1)
  block <- paste(CODE[at:(at + 3)], collapse = " ")
  expect_true(grepl("check.names = FALSE", block, fixed = TRUE))
})

test_that("a curated exclusion is not reinstated by the vocabulary route", {
  # {acc}_groups.csv is the curator's contrast. On the supplementary-count
  # path the vocabulary route re-derives labels from column names, and a
  # column name cannot express an exclusion: GSE217932's curated 4 v 4 became
  # 8 v 4 again because the THP-treated arm is named like a case, and the log
  # said so plainly -- "Manual annotation: 4 case, 4 control, 4 excluded"
  # immediately followed by "Column vocab: 8 case, 4 control". The title and
  # metadata routes carry pheno_used$group, which is the curated label, so
  # when a curator has excluded samples the vocabulary route must be skipped.
  at <- grep("col_groups <- classify_cols_by_vocab(colnames(mat_num))",
             CODE, fixed = TRUE)
  expect_length(at, 1)
  block <- paste(CODE[max(1, at - 8):(at + 2)], collapse = " ")
  expect_true(grepl("isTRUE(curated)", block, fixed = TRUE))

  # And the flag has to reach the function, or the guard above is unreachable.
  sig <- grep("run_deseq2_from_matrix <- function", CODE, fixed = TRUE)
  expect_length(sig, 1)
  expect_true(grepl("curated", paste(CODE[sig:(sig + 1)], collapse = " ")))
  calls <- grep("return(run_deseq2_from_matrix(", CODE, fixed = TRUE)
  expect_gt(length(calls), 0)
  for (i in calls) {
    expect_true(grepl("curated = curated_exclusions",
                      paste(CODE[i:(i + 1)], collapse = " "), fixed = TRUE))
  }
})

test_that("an ambiguous title match is disambiguated at a word boundary", {
  # Containment alone cannot separate labels where one contains another, and
  # the plainest case in this corpus is "male" inside "female": every Male
  # column of GSE197233 matched two titles, its own and the Female one, and
  # was dropped as ambiguous -- halving the study without a word of complaint.
  at <- grep("matched_cols <- vapply(cols_lo,", CODE, fixed = TRUE)
  expect_length(at, 1)
  block <- paste(CODE[at:(at + 6)], collapse = " ")
  expect_true(grepl("at_boundary", block, fixed = TRUE))
  expect_true(any(grepl("at_boundary <- function", CODE, fixed = TRUE)))
})

test_that("a derived unit takes its phenotype from the parent accession", {
  # GEO knows GSE180627, not GSE180627_PAG. Everything else stays keyed on the
  # derived id; only getGEO sees the parent.
  at <- grep('geo_accession <- sub("_.*$", "", accession)', CODE, fixed = TRUE)
  expect_length(at, 1)
  expect_true(any(grepl("getGEO(geo_accession", CODE, fixed = TRUE)))
})

test_that("the title match takes an exact hit before any containment hit", {
  # A word boundary does not separate "mn_1" from "mn_10": the column matches
  # both titles at position 1. Only exactness does. GSE250152 lost 18 of its
  # 33 samples to this and GSE186505 two of twenty, and both looked like
  # studies whose metadata simply did not line up.
  at <- grep("exact <- which(titles == nm)", CODE, fixed = TRUE)
  expect_length(at, 1)
  # and it must return before the containment search, not after it
  tail <- paste(CODE[at:(at + 2)], collapse = " ")
  expect_true(grepl("return(pheno_used$group[exact])", tail, fixed = TRUE))
  containment <- grep("idx <- which(sapply(titles,", CODE, fixed = TRUE)
  expect_true(all(containment > at))
})
