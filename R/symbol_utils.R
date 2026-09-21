#' Gene symbol canonicalisation.
#'
#' Two failure modes fragment a gene across studies, and both are invisible
#' downstream because every join in the pipeline is a case-sensitive string
#' match on the symbol:
#'
#'   1. Excel date corruption. Spreadsheet software parses MARCH1 and SEPT1 as
#'      dates and writes them back as "1-Mar" and "1-Sep". The damage is done
#'      in the submitted supplementary files, so it arrives with the data.
#'   2. Superseded nomenclature. HGNC renamed MARCH1-11 to MARCHF1-11 and
#'      SEPT1-14 to SEPTIN1-14, so studies annotated at different times use
#'      different symbols for the same gene.
#'
#' Left alone, SEPTIN1 pools three times over: as Septin1, Sept1 and 1-Sep.
#'
#' Scope is deliberately narrow. Only the MARCH/MARCHF and SEPT/SEPTIN
#' families, plus SEP15, are touched: these are the families Excel corrupts
#' and their renaming is unambiguous. Notably 1-Dec is *not* repaired --- DEC1
#' is an alias for both BHLHE40 and DELEC1, so it cannot be resolved from the
#' symbol alone, and no observed study carries it.

# MARCH1-11 are the complete MARCHF family.
.MARCHF_N <- 1:11

# Septins are SEPT1-SEPT12 and SEPT14. There is no SEPT13, and 15 is not a
# septin at all: "15-Sep" is SEP15, the 15 kDa selenoprotein, now SELENOF.
# Mapping it to a septin would invent a gene that does not exist.
.SEPTIN_N <- c(1:12, 14)
.SELENOF_N <- 15


#' Detect whether a table spells gene symbols in upper or title case.
#'
#' Human tables use HGNC all-caps (TNF); rodent tables use title case (Tnf).
#' Bare numeric probe ids carry no case information and are ignored, or a
#' rodent array table would be decided by its probe ids rather than its genes.
#'
#' @param ids Character vector of feature ids.
#' @return "upper" or "title". Ties and absent evidence give "upper", the
#'   HGNC convention.
symbol_case_style <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  has_alpha <- grepl("[A-Za-z]", ids)
  ids <- ids[has_alpha]
  if (length(ids) == 0) return("upper")

  n_upper <- sum(grepl("^[A-Z][A-Z0-9._-]*$", ids))
  n_title <- sum(grepl("^[A-Z][a-z]", ids))
  if (n_title > n_upper) "title" else "upper"
}


#' Apply a case style to a canonical symbol given in upper case.
.apply_style <- function(symbol_upper, style) {
  if (identical(style, "title")) {
    paste0(substr(symbol_upper, 1, 1),
           tolower(substr(symbol_upper, 2, nchar(symbol_upper))))
  } else {
    symbol_upper
  }
}


#' Resolve one id to its canonical symbol, or NA when it is not in scope.
.canonical_or_na <- function(id, style) {
  if (is.na(id) || !nzchar(id)) return(NA_character_)

  # Excel dates, in either order: "1-Mar" / "Mar-1".
  m <- regmatches(id, regexec("^([0-9]{1,2})-(Mar|Sep)$", id, ignore.case = TRUE))[[1]]
  if (length(m) == 0) {
    m <- regmatches(id, regexec("^(Mar|Sep)-([0-9]{1,2})$", id, ignore.case = TRUE))[[1]]
    if (length(m) == 3) m <- c(m[1], m[3], m[2])
  }
  if (length(m) == 3) {
    n <- as.integer(m[2])
    month <- tolower(m[3])
    if (month == "mar" && n %in% .MARCHF_N)
      return(.apply_style(paste0("MARCHF", n), style))
    if (month == "sep" && n %in% .SEPTIN_N)
      return(.apply_style(paste0("SEPTIN", n), style))
    if (month == "sep" && n %in% .SELENOF_N)
      return(.apply_style("SELENOF", style))
    return(NA_character_)   # e.g. 13-Sep, 12-Mar: not attributable
  }

  # Current and legacy family names, any casing.
  m <- regmatches(id, regexec("^(MARCHF|MARCH|SEPTIN|SEPT|SEP)([0-9]{1,2})$",
                              id, ignore.case = TRUE))[[1]]
  if (length(m) == 3) {
    prefix <- toupper(m[2])
    n <- as.integer(m[3])
    if (prefix %in% c("MARCH", "MARCHF") && n %in% .MARCHF_N)
      return(.apply_style(paste0("MARCHF", n), style))
    if (prefix %in% c("SEPT", "SEPTIN") && n %in% .SEPTIN_N)
      return(.apply_style(paste0("SEPTIN", n), style))
    if (prefix == "SEP" && n %in% .SELENOF_N)
      return(.apply_style("SELENOF", style))
    return(NA_character_)
  }

  if (identical(toupper(id), "SELENOF")) return(.apply_style("SELENOF", style))

  NA_character_
}


#' Keep one row per feature_id, the one with the smallest p-value.
#'
#' Base-R equivalent of the group_by/slice_min rule aggregate_to_gene uses.
#' Original row order is preserved so a repaired table diffs minimally against
#' the one it replaces.
keep_min_pval_per_feature <- function(df) {
  if (!"pval" %in% colnames(df)) return(df[!duplicated(df$feature_id), , drop = FALSE])
  o <- order(df$feature_id, df$pval, na.last = TRUE)
  first <- o[!duplicated(df$feature_id[o])]
  df[sort(first), , drop = FALSE]
}


#' Repair a study's symbol column and re-aggregate any genes the repair merges.
#'
#' Shared by 05b and scripts/migrate_canonicalise_symbols.R so the two cannot
#' drift. Re-aggregation is the part that is easy to forget: a study carrying
#' both 1-Sep and Sept1 collapses to two Septin1 rows, which a meta-analysis
#' would read as two independent observations of one gene from one study.
#'
#' @param df Data frame with a feature_id column.
#' @param accession Study id, used only in messages.
#' @param verbose Whether to report what changed.
#' @return List of `df` (repaired) and `changes` (data frame of from/to pairs).
canonicalize_study_symbols <- function(df, accession = "", verbose = TRUE) {
  empty <- data.frame(study_id = character(0), from = character(0),
                      to = character(0), stringsAsFactors = FALSE)
  if (!"feature_id" %in% colnames(df) || nrow(df) == 0) {
    return(list(df = df, changes = empty))
  }

  original <- as.character(df$feature_id)
  style    <- symbol_case_style(original)
  repaired <- canonicalize_gene_symbols(original, style)

  changed <- which(!is.na(repaired) & !is.na(original) & repaired != original)
  if (length(changed) == 0) return(list(df = df, changes = empty))

  changes <- unique(data.frame(
    study_id = accession, from = original[changed], to = repaired[changed],
    stringsAsFactors = FALSE
  ))

  if (verbose) {
    ex <- utils::head(paste0(changes$from, "->", changes$to), 4)
    message(sprintf("  %s: canonicalised %d row(s) [%s style]: %s%s",
                    accession, length(changed), style, paste(ex, collapse = ", "),
                    if (nrow(changes) > length(ex)) ", ..." else ""))
  }

  df$feature_id <- repaired
  if (!"original_feature_id" %in% colnames(df)) df$original_feature_id <- original

  dup <- sum(duplicated(df$feature_id))
  if (dup > 0) {
    if (verbose) {
      message(sprintf("  %s: %d row(s) merged (min pval kept)", accession, dup))
    }
    df <- keep_min_pval_per_feature(df)
  }
  list(df = df, changes = changes)
}


#' Canonicalise gene symbols, repairing Excel dates and legacy family names.
#'
#' Ids outside the MARCH/SEPT families are returned unchanged, including NA and
#' empty strings. Two different ids can map to the same symbol (a study holding
#' both 1-Sep and Sept1), so callers must re-aggregate afterwards or the study
#' will contribute one gene to a meta-analysis twice.
#'
#' @param ids Character vector of feature ids.
#' @param style "upper" or "title"; defaults to the majority style in `ids`.
#' @return Character vector the same length as `ids`.
canonicalize_gene_symbols <- function(ids, style = NULL) {
  if (length(ids) == 0) return(ids)
  if (is.null(style)) style <- symbol_case_style(ids)

  ids <- as.character(ids)
  # Only a handful of distinct ids are ever in scope, so resolve over the
  # unique set rather than per row.
  uniq <- unique(ids)
  mapped <- vapply(uniq, .canonical_or_na, character(1), style = style,
                   USE.NAMES = FALSE)
  names(mapped) <- uniq

  out <- ids
  hit <- !is.na(mapped[ids])
  out[hit] <- mapped[ids][hit]
  out
}
