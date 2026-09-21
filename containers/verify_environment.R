#!/usr/bin/env Rscript
#' Prove the container closed the reproducibility gaps, rather than asserting it.
#'
#' Two failures this checks for are silent by nature, which is the whole reason
#' the check exists:
#'
#'   - A missing annotation package does not make 05b fail. It falls through to
#'     "keeping native IDs", nine studies degrade from gene symbols to probe
#'     ids, and every downstream number moves for a reason nobody sees.
#'   - `_dependencies.R` exists so renv::snapshot() can see packages 05b selects
#'     from a lookup table. It listed 5 of the 11 for months. A lockfile built
#'     from an incomplete shim installs an incomplete environment, and the first
#'     symptom is wrong output rather than a crash.
#'
#' So this compares GPL_TO_ANNO_PKG against _dependencies.R against what is
#' actually loadable, and fails if the three disagree.
#'
#' Usage: Rscript containers/verify_environment.R

# Walk up from the working directory to the repository root, so the script runs
# the same from /work in the container and from a checkout on a workstation.
# `%||%` is not used here: base R gained it only in 4.4, and this script must
# run on the older interpreter it is diagnosing.
find_repo <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    if (dir.exists(file.path(dir, "pipeline")) &&
        file.exists(file.path(dir, "renv.lock"))) return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    dir <- parent
  }
  stop("could not find the repository root (no pipeline/ + renv.lock above ", start, ")")
}
repo <- find_repo()

ok <- TRUE
fail <- function(...) { message("FAIL  ", ...); ok <<- FALSE }
pass <- function(...) message("ok    ", ...)

# ---------------------------------------------------------------------------
# 1. R and Bioconductor versions match what renv.lock pins
# ---------------------------------------------------------------------------
lock <- jsonlite::fromJSON(file.path(repo, "renv.lock"))
want_r <- lock$R$Version
got_r <- as.character(getRversion())
if (substr(got_r, 1, 3) == substr(want_r, 1, 3)) {
  pass(sprintf("R %s (lockfile pins %s)", got_r, want_r))
} else {
  fail(sprintf("R %s but lockfile pins %s: renv::restore() will not resolve",
               got_r, want_r))
}

want_bioc <- lock$Bioconductor$Version
got_bioc <- as.character(BiocManager::version())
if (identical(substr(got_bioc, 1, 4), substr(want_bioc, 1, 4))) {
  pass(sprintf("Bioconductor %s (lockfile pins %s)", got_bioc, want_bioc))
} else {
  fail(sprintf("Bioconductor %s but lockfile pins %s", got_bioc, want_bioc))
}

# ---------------------------------------------------------------------------
# 2. Every annotation package 05b can select is installed
# ---------------------------------------------------------------------------
# Parsed from the source rather than restated, so this cannot drift from the
# table it is checking.
src <- readLines(file.path(repo, "pipeline", "05b_normalize_ids.R"), warn = FALSE)
start <- grep("^GPL_TO_ANNO_PKG <- list\\(", src)
if (length(start) != 1) {
  fail("could not locate GPL_TO_ANNO_PKG in pipeline/05b_normalize_ids.R")
  anno <- character(0)
} else {
  end <- start + which(grepl("^\\)", src[(start + 1):length(src)]))[1]
  block <- src[start:end]
  anno <- unique(gsub('.*"([^"]+\\.db)".*', "\\1",
                      grep('"[^"]+\\.db"', block, value = TRUE)))
  pass(sprintf("GPL_TO_ANNO_PKG names %d distinct annotation packages", length(anno)))
}

missing <- anno[!vapply(anno, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  fail(sprintf("%d annotation package(s) absent: %s", length(missing),
               paste(missing, collapse = ", ")))
  fail("  05b would not error -- it would silently keep probe ids as feature ids")
} else if (length(anno)) {
  pass(sprintf("all %d annotation packages loadable", length(anno)))
}

# ---------------------------------------------------------------------------
# 3. _dependencies.R covers the same set, so future snapshots stay complete
# ---------------------------------------------------------------------------
dep_src <- readLines(file.path(repo, "_dependencies.R"), warn = FALSE)
declared <- gsub(".*library\\(([^)]+)\\).*", "\\1",
                 grep("^library\\(", dep_src, value = TRUE))
undeclared <- setdiff(anno, declared)
if (length(undeclared)) {
  fail(sprintf("in GPL_TO_ANNO_PKG but not _dependencies.R: %s",
               paste(undeclared, collapse = ", ")))
  fail("  renv::snapshot() cannot see these, so the lockfile will omit them")
} else if (length(anno)) {
  pass("_dependencies.R covers every entry of GPL_TO_ANNO_PKG")
}

# ---------------------------------------------------------------------------
# 4. The analysis packages 05/05b/06 need
# ---------------------------------------------------------------------------
core <- c("DESeq2", "edgeR", "limma", "biomaRt", "metafor", "dplyr",
          "data.table", "jsonlite", "yaml", "httr", "testthat",
          "org.Mm.eg.db", "org.Rn.eg.db")
absent <- core[!vapply(core, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  fail(sprintf("core packages absent: %s", paste(absent, collapse = ", ")))
} else {
  pass(sprintf("all %d core analysis packages loadable", length(core)))
}

# ---------------------------------------------------------------------------
# 5. MAGMA
# ---------------------------------------------------------------------------
magma <- Sys.which("magma")
if (nzchar(magma)) {
  ver <- tryCatch(system2("magma", "--version", stdout = TRUE, stderr = TRUE),
                  error = function(e) "")
  pass(sprintf("magma on PATH: %s", paste(head(ver, 1), collapse = " ")))
} else {
  fail("magma not on PATH: the genomics arm cannot be re-run")
}

ref <- file.path(repo, "tools", "magma", "g1000_eur.bed")
if (file.exists(ref)) {
  pass("1000G EUR reference present")
} else {
  message("note  1000G EUR reference absent -- mount it or run ",
          "containers/fetch_reference_data.sh (not a failure: it is a 4 GB ",
          "volume, deliberately outside the image)")
}

message("")
if (ok) {
  message("environment verified: 05b can map probes to symbols, and the ",
          "genomics arm has its tooling")
  quit(status = 0)
}
message("environment INCOMPLETE -- see FAIL lines above")
quit(status = 1)
