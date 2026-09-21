#' Shared utilities for per-study differential analysis (step 05).
#'
#' Covers: harmonized manifest loading, sample-group labelling,
#' effect-size standardization, and per-study output writing.

library(jsonlite)
library(yaml)
library(dplyr)


# ---------------------------------------------------------------------------
# Manifest I/O
# ---------------------------------------------------------------------------

#' Load harmonized dataset records for one modality.
#'
#' @param interim_dir Path to data/interim/.
#' @param modality    Character string: one of the five modalities.
#' @return data.frame with one row per dataset.
load_harmonized <- function(interim_dir, modality) {
  path <- file.path(interim_dir, modality, "harmonized.jsonl")
  if (!file.exists(path)) {
    message(sprintf("[%s] No harmonized manifest at %s", modality, path))
    return(data.frame())
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  # Parse each JSON record; collapse list-typed fields to semicolon-separated strings
  records <- lapply(lines, function(line) {
    rec <- fromJSON(line, simplifyVector = TRUE)
    for (nm in names(rec)) {
      if (is.list(rec[[nm]]) || length(rec[[nm]]) > 1) {
        rec[[nm]] <- paste(unlist(rec[[nm]]), collapse = ";")
      }
    }
    as.data.frame(rec, stringsAsFactors = FALSE)
  })
  df <- bind_rows(records)
  message(sprintf("[%s] Loaded %d harmonized records", modality, nrow(df)))
  df
}


#' Screening records that are keyed by dataset accession, in preference order.
#'
#' NOT literature/prisma/screening.csv: that is the step 02 literature sheet,
#' keyed by PubMed id. Matching its `uid` column against GEO/PRIDE/
#' MetaboLights accessions compares two different identifier spaces and
#' narrows the corpus to whatever accessions happen to have been typed into
#' the sheet — 13 of 283 transcriptomic records, 0 of every other modality.
#' See plan/2026-08-28-step05-screening-mismatch.md.
SCREENING_RECORD_SUFFIXES <- c("_candidates.csv", "_manual_review.csv")


#' Leading alphabetic run of a dataset accession, e.g. "GSE241361_DRG" -> "GSE".
#'
#' Named for accessions specifically: `id_space` is already used across this
#' codebase as a variable and column name for feature-identifier spaces
#' ("HGNC", "ChEBI"), which is a different thing.
#'
#' PubMed ids are all digits and so map to "", which is what makes them
#' distinguishable from every repository accession space.
#'
#' @param x Character vector of identifiers.
#' @return Character vector of identifier-space prefixes.
accession_id_space <- function(x) {
  sub("[^A-Za-z].*$", "", as.character(x))
}


#' Load one modality's screening record and return its included accessions.
#'
#' Reads the first of literature/prisma/{modality}_candidates.csv or
#' literature/prisma/{modality}_manual_review.csv that exists, and returns the
#' `accession` values whose `verdict` is "include". Both files are maintained
#' in dataset-accession space, unlike the literature sheet.
#'
#' @param prisma_dir Path to literature/prisma/.
#' @param modality   Character string: one of the five modalities.
#' @return Character vector of included accessions, or NULL when the modality
#'   has no screening record (caller should then process every harmonized
#'   record).
load_included_accessions <- function(prisma_dir, modality) {
  paths <- file.path(prisma_dir, paste0(modality, SCREENING_RECORD_SUFFIXES))
  path  <- paths[file.exists(paths)][1]

  if (is.na(path)) {
    message(sprintf(
      "[%s] No screening record (%s) — including all harmonized records.",
      modality, paste(basename(paths), collapse = " or ")))
    return(NULL)
  }

  df <- read.csv(path, stringsAsFactors = FALSE)
  missing_cols <- setdiff(c("accession", "verdict"), colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("%s lacks required column(s): %s",
                 path, paste(missing_cols, collapse = ", ")))
  }

  # %in%, not ==: a missing verdict must not put an NA into the accession list.
  included <- df$accession[df$verdict %in% "include"]
  message(sprintf("[%s] Screening: %d of %d records included (%s)",
                  modality, length(included), nrow(df), basename(path)))
  included
}


#' Which harmonized accessions a screening record admits.
#'
#' Exact match, or -- for a derived unit -- the verdict of the parent
#' accession it was split from. Only parents are ever screened: a derived unit
#' (GSE180627_S1, GSE241361_DRG) is created by the analysis, not retrieved, so
#' it appears in no screening sheet by construction. A plain `%in%` therefore
#' drops every one of them, and does it the way this filter has failed before
#' -- by producing a smaller corpus of exactly the right shape. On 2026-08-31
#' that was six of 83 records, the whole of the GSE180627 and GSE197233
#' regional splits, with `05` reporting "77 studies to process" and no warning.
#'
#' Suffixes are stripped one `_<token>` at a time so a unit derived twice still
#' resolves. `_<token>` is not part of any accession space this project reads
#' (GSE, PXD, MTBLS, ST), so nothing but a derived id can match.
#'
#' @param accessions  Character vector from the harmonized manifest.
#' @param included_ids Character vector from load_included_accessions().
#' @return Logical vector the length of `accessions`.
accession_screened_in <- function(accessions, included_ids) {
  accessions <- as.character(accessions)
  ok <- accessions %in% included_ids
  candidate <- accessions
  repeat {
    miss <- which(!ok)
    if (length(miss) == 0) break
    parent <- sub("_[^_]+$", "", candidate[miss])
    if (all(parent == candidate[miss])) break
    candidate[miss] <- parent
    ok[miss] <- parent %in% included_ids
  }
  ok
}


#' Stop unless the included accessions can actually select harmonized records.
#'
#' The defect this guards against narrowed the corpus without erroring: a
#' non-empty list in the wrong identifier space produces the same output shape
#' as a correct one, just with almost every study missing. Two conditions are
#' checked, both of which the old screening.csv path violated:
#'
#' 1. Most included ids must share an identifier space with the manifest. A
#'    prefix-intersection test alone is not enough — 13 of the 47 ids in the
#'    literature sheet were GEO accessions, so the spaces did overlap while
#'    the file was still the wrong one.
#' 2. At least one included accession must name a harmonized record.
#'
#' @param included_ids Character vector from load_included_accessions(), or NULL.
#' @param accessions   Character vector of harmonized accessions.
#' @param modality     Character string, for the error message.
#' @param min_shared   Minimum fraction of included ids that must lie in the
#'   manifest's identifier space.
#' @return Invisibly TRUE; stops otherwise.
assert_screening_identifier_space <- function(included_ids, accessions,
                                              modality, min_shared = 0.5) {
  if (is.null(included_ids) || length(included_ids) == 0) return(invisible(TRUE))

  manifest_spaces <- unique(accession_id_space(accessions))
  shared <- mean(accession_id_space(included_ids) %in% manifest_spaces)
  if (shared < min_shared) {
    stop(sprintf(
      paste0("[%s] screening record is in the wrong identifier space: only ",
             "%d of %d included ids (%.0f%%) look like harmonized accessions. ",
             "Included ids look like: %s. Harmonized accessions look like: %s."),
      modality, sum(accession_id_space(included_ids) %in% manifest_spaces),
      length(included_ids), 100 * shared,
      paste(utils::head(unique(included_ids), 3), collapse = ", "),
      paste(utils::head(unique(accessions), 3), collapse = ", ")))
  }

  if (length(intersect(included_ids, accessions)) == 0) {
    stop(sprintf(
      paste0("[%s] no harmonized record matches any of the %d included ",
             "accessions — refusing to process zero studies."),
      modality, length(included_ids)))
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Sample grouping
# ---------------------------------------------------------------------------

#' Assign case/control labels to a GEO phenotype data.frame.
#'
#' Scans `characteristics_ch1` (or `title`) for case_labels / control_labels
#' defined in the analysis config.
#'
#' @param pheno       data.frame from GEOquery pData().
#' @param case_labels Character vector of case substrings.
#' @param ctrl_labels Character vector of control substrings.
#' @return pheno with an added `group` column ("case" | "control" | NA).
assign_groups <- function(pheno, case_labels, ctrl_labels) {
  # If pheno already has a pre-assigned group column, use it directly
  if ("group" %in% colnames(pheno) &&
      all(pheno$group %in% c("case", "control", NA))) {
    n_case <- sum(pheno$group == "case", na.rm = TRUE)
    n_ctrl <- sum(pheno$group == "control", na.rm = TRUE)
    message(sprintf("  Group assignment: %d case, %d control, %d unassigned",
                    n_case, n_ctrl, sum(is.na(pheno$group))))
    return(pheno)
  }

  search_cols <- intersect(
    c("characteristics_ch1", "characteristics_ch1.1", "title", "source_name_ch1"),
    colnames(pheno)
  )
  if (length(search_cols) == 0) {
    pheno$group <- NA_character_
    message("  Group assignment: 0 case, 0 control, no searchable columns")
    return(pheno)
  }
  combined <- apply(pheno[, search_cols, drop = FALSE], 1, function(row) {
    tolower(paste(row, collapse = " "))
  })

  group <- rep(NA_character_, nrow(pheno))

  # A control arm routinely names the disease it does not have -- "no
  # endometriosis", "non-endometriosis", "disease-free". Matching is fixed
  # substring and first-match-wins, so the case label hits that value first and
  # the control becomes a case, leaving a study with two case arms and no
  # split. Negation is therefore resolved before either loop, the same way
  # pipeline/auto_annotate_geo.py does it. Built by concatenation rather than
  # regex so the labels need no escaping.
  negated <- rep(FALSE, length(combined))
  neg_terms <- unlist(lapply(
    tolower(c(case_labels, "disease", "pain", "lesion", "symptoms")),
    function(l) c(paste0("non", l), paste0("non-", l), paste0("non_", l),
                  paste0("non ", l), paste0("no ", l), paste0("without ", l),
                  paste0(l, "-free"), paste0(l, " free"))))
  for (term in neg_terms) {
    negated <- negated | grepl(term, combined, fixed = TRUE)
  }
  group[negated] <- "control"

  for (lbl in tolower(case_labels)) {
    group[grepl(lbl, combined, fixed = TRUE) & is.na(group)] <- "case"
  }
  for (lbl in tolower(ctrl_labels)) {
    group[grepl(lbl, combined, fixed = TRUE) & is.na(group)] <- "control"
  }

  pheno$group <- group
  n_case <- sum(group == "case", na.rm = TRUE)
  n_ctrl <- sum(group == "control", na.rm = TRUE)
  message(sprintf("  Group assignment: %d case, %d control, %d unassigned",
                  n_case, n_ctrl, sum(is.na(group))))
  pheno
}


# ---------------------------------------------------------------------------
# Effect-size output format
# ---------------------------------------------------------------------------

#' Build a standardized per-study effect-size data.frame.
#'
#' All downstream steps expect exactly these columns.
#'
#' @param feature_id    Character vector of gene/metabolite IDs.
#' @param effect_size   Numeric effect size (log2FC, beta, Hedges g, SMD).
#' @param se            Standard error of effect_size.
#' @param pval          Raw p-value.
#' @param padj          Adjusted p-value.
#' @param n_case        Number of case samples.
#' @param n_control     Number of control samples.
#' @param study_id      GEO/PRIDE/MetaboLights accession.
#' @param modality      Modality string.
#' @param id_space      Identifier namespace (HGNC, UniProt, ChEBI, …).
#' @param effect_unit   Unit string (log2FC, beta, SMD, OR).
#' @return data.frame with standardized columns.
build_effect_df <- function(feature_id, effect_size, se, pval, padj,
                            n_case, n_control, study_id, modality,
                            id_space, effect_unit) {
  data.frame(
    feature_id   = as.character(feature_id),
    effect_size  = as.numeric(effect_size),
    se           = as.numeric(se),
    pval         = as.numeric(pval),
    padj         = as.numeric(padj),
    n_case       = as.integer(n_case),
    n_control    = as.integer(n_control),
    study_id     = as.character(study_id),
    modality     = as.character(modality),
    id_space     = as.character(id_space),
    effect_unit  = as.character(effect_unit),
    stringsAsFactors = FALSE
  )
}


#' Write per-study effect data.frame to results/per_study/{modality}/.
#'
#' @param df         data.frame from build_effect_df().
#' @param out_dir    Base per_study output directory.
#' @param modality   Modality string.
#' @param accession  Study accession (used as filename prefix).
save_effects <- function(df, out_dir, modality, accession) {
  mod_dir <- file.path(out_dir, modality)
  dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(mod_dir, paste0(accession, "_effects.csv"))
  write.csv(df, out_path, row.names = FALSE)
  message(sprintf("  Saved %d effect-size rows → %s", nrow(df), out_path))
  invisible(out_path)
}
