#' Per-study differential analysis for transcriptomics (RNA-seq and microarray).
#'
#' RNA-seq:    DESeq2 with Wald test, using the UNSHRUNKEN results() standard
#'             error. Do not substitute lfcShrink(type = "apeglm"): its posterior
#'             SD (~0.001) is two orders of magnitude below the true SE
#'             (0.10-0.25) and would corrupt every pooled estimate.
#' Microarray: limma with eBayes; SE = |logFC / t| (see R/se_utils.R).
#'
#' Requires: GEOquery, DESeq2, limma, edgeR, AnnotationDbi, org.Hs.eg.db

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(dplyr)
})

# Directory of *this* file, for sourcing its siblings.
#
# These files live in R/ but are sourced by scripts in pipeline/, and the
# previous idiom derived the directory from commandArgs("--file="), which names
# the *invoking* script. Called as `Rscript pipeline/05_per_study_da.R` that
# resolved to pipeline/se_utils.R, which does not exist. It went unnoticed
# because each source() is guarded by an `exists()` check, so it only fires
# when nothing has already defined the function and the failure depends on
# source order.
#
# When a file is source()d, R records its path as `ofile` on the sourcing
# frame; walking the stack finds it. The commandArgs fallback still applies
# when this file is itself the script being run.
.sibling_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("--file=", "", fa[1]))))
  getwd()
}

if (!exists("assign_groups")) source(file.path(.sibling_dir(), "da_utils.R"))
if (!exists("se_from_pvalue")) source(file.path(.sibling_dir(), "se_utils.R"))
if (!exists("run_cuffdiff_output")) source(file.path(.sibling_dir(), "cuffdiff_da.R"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

#' Run differential analysis for one GEO transcriptomics study.
#'
#' @param accession  GEO series accession (e.g. "GSE12345").
#' @param platform   Canonical platform string ("RNA-seq" or "microarray").
#' @param cfg        Named list from conf/analysis/default.yaml.
#' @param cache_dir  Directory for GEO data cache (passed to getGEO).
#' @return data.frame of effect sizes, or NULL on failure.
run_transcriptomics_da <- function(accession, platform, cfg, cache_dir = tempdir()) {
  message(sprintf("\n=== %s [%s] ===", accession, platform))

  # A derived unit (GSE180627_PAG, GSE241361_DRG) is one region of a parent
  # accession, so GEO knows only the parent. Everything else -- the curated
  # groups file, the supplementary directory, the output name -- stays keyed
  # on the derived id; only the phenotype comes from the parent, which is
  # exactly right, because the derived groups.csv is what restricts it.
  geo_accession <- sub("_.*$", "", accession)
  if (geo_accession != accession) {
    message("  Derived unit: taking the phenotype from ", geo_accession)
  }
  geo_data <- tryCatch(
    getGEO(geo_accession, GSEMatrix = TRUE, destdir = cache_dir, getGPL = FALSE),
    error = function(e) { message("  GEO download failed: ", e$message); NULL }
  )
  if (is.null(geo_data)) return(NULL)

  eset <- geo_data[[1]]
  pheno <- pData(eset)

  # Per-study manual annotation overrides automatic keyword matching
  manual_path <- file.path(cache_dir, paste0(accession, "_groups.csv"))
  if (file.exists(manual_path)) {
    message("  Loading manual sample annotation: ", basename(manual_path))
    ann <- read.csv(manual_path, stringsAsFactors = FALSE)
    # ann has columns: sample_id, group ("case" | "control" | "")
    pheno$group <- ann$group[match(rownames(pheno), ann$sample_id)]
    pheno$group[pheno$group == ""] <- NA
    curated_exclusions <- any(is.na(pheno$group)) && any(!is.na(pheno$group))
    n_pre <- sum(!is.na(pheno$group))
    message(sprintf("  Manual annotation: %d case, %d control, %d excluded",
                    sum(pheno$group == "case", na.rm = TRUE),
                    sum(pheno$group == "control", na.rm = TRUE),
                    sum(is.na(pheno$group))))
  } else {
    curated_exclusions <- FALSE
    pheno <- assign_groups(pheno,
                           case_labels = cfg$da$case_labels,
                           ctrl_labels  = cfg$da$control_labels)
  }

  # Where the expression matrix comes from decides whether the series-matrix
  # phenotype may gate this study at all, so resolve it first. A supplementary
  # count file carries its own sample columns, and run_deseq2_from_matrix()
  # groups them with classify_cols_by_vocab(), consulting this phenotype only
  # as a fallback. Gating on it here dropped ten studies whose grouping is
  # perfectly well defined one step later -- and because save_effects() then
  # never fires, their stale tables silently survived. The CEL branch does use
  # `usable`, `n_case` and `n_ctrl`, so it stays gated.
  supp_dir  <- file.path(cache_dir, accession)
  cel_files <- list.files(supp_dir, pattern = "\\.CEL\\.gz$",
                          ignore.case = TRUE, full.names = TRUE)
  suppl_file <- find_suppl_count_file(supp_dir)

  usable <- pheno[!is.na(pheno$group), , drop = FALSE]
  n_case <- sum(usable$group == "case")
  n_ctrl <- sum(usable$group == "control")

  if (is.null(suppl_file) &&
      (n_case < cfg$da$min_samples_per_group ||
       n_ctrl < cfg$da$min_samples_per_group)) {
    message(sprintf("  Skipping: insufficient group sizes (case=%d, ctrl=%d)", n_case, n_ctrl))
    return(NULL)
  }

  if (!is.null(suppl_file)) {
    ext <- tolower(suppl_file)
    message("  Using supplementary file: ", basename(suppl_file))
    if (grepl("\\.xlsx?$", ext)) {
      return(run_deseq2_from_matrix(suppl_file, usable, accession, cfg,
                                   curated = curated_exclusions))
    } else {
      return(run_deseq2_from_matrix(suppl_file, usable, accession, cfg,
                                   curated = curated_exclusions))
    }
  } else if (length(cel_files) > 0) {
    message("  Using CEL files (", length(cel_files), " files) for microarray normalization")
    return(run_limma_from_cel(cel_files, usable, n_case, n_ctrl, accession, cfg))
  } else if (nrow(exprs(eset)) == 0) {
    message("  No expression data in series matrix or supplementary files — skipping")
    return(NULL)
  } else if (platform == "RNA-seq") {
    run_deseq2(eset, usable, n_case, n_ctrl, accession, cfg)
  } else {
    run_limma_array(eset, usable, n_case, n_ctrl, accession, cfg)
  }
}


# ---------------------------------------------------------------------------
# Supplementary file discovery
# ---------------------------------------------------------------------------

COUNT_PATTERNS_R <- paste0(
  "count|raw.count|gene.count|expression.matrix|gene.exp|read.count",
  "|feature.count|htseq|deseq|rna.seq|counts.matrix",
  "|count_matrix|raw_counts|gene_counts|read_counts",
  "|tpm|rpkm|fpkm|cpm|normalized.count",
  "|expression.gene|gene.expression|read.per.feature|reads.per",
  "|matrix\\.mtx|sparse.matrix|allsamples.*count|combined.*count",
  "|\\.count\\.|\\.counts\\.|expression_data|exprs|expr_mat"
)

SKIP_PATTERNS_R <- "\\.(bam|bai|cram|fastq|fq|bigwig|bw|bed|gtf|gff|vcf|sra)(\\.gz)?$"

#' Find the best supplementary count/expression matrix file in a directory.
#' Returns the file path or NULL.
find_suppl_count_file <- function(dir_path) {
  if (!dir.exists(dir_path)) return(NULL)
  files <- list.files(dir_path, full.names = TRUE)
  # Filter out skip patterns and metadata/annotation sidecar files
  files <- files[!grepl(SKIP_PATTERNS_R, files, ignore.case = TRUE)]
  files <- files[!grepl("_metadata\\.|_meta\\.|_annotation\\.|_coldata\\.",
                         basename(files), ignore.case = TRUE)]

  # Priority tiers (earlier = preferred):
  # 1. Raw integer count files
  raw_pat  <- "raw[._-]?(umi[._-]?)?count|htseq|read[._-]?count|umi[._-]?count|_counts\\.(txt|tsv|csv|xlsx?)(|\\.gz)$"
  # 2. Cuffdiff / pre-computed DE output (treated as effect size input)
  cuffdiff_pat <- "gene_expression\\.tsv|cuffdiff|diff_exp|differential"
  # 3. Normalized expression (TPM/FPKM)
  norm_pat <- "tpm|fpkm|rpkm|cpm|normalized"
  xlsx_pat <- "\\.xlsx?$"

  raw_files      <- files[grepl(raw_pat,      basename(files), ignore.case = TRUE)]
  cuffdiff_files <- files[grepl(cuffdiff_pat, basename(files), ignore.case = TRUE)]
  norm_files     <- files[grepl(norm_pat,     basename(files), ignore.case = TRUE)]
  xlsx_files     <- files[grepl(xlsx_pat,     basename(files), ignore.case = TRUE) &
                           !grepl(norm_pat,   basename(files), ignore.case = TRUE)]
  any_files      <- files[grepl(COUNT_PATTERNS_R, basename(files), ignore.case = TRUE)]

  # Prefer: raw > cuffdiff > xlsx > norm > any. Deliberately no preference by
  # filename within the raw tier: 05d's merged matrix looked like the safer
  # choice, but for GSE135080 its merge dropped the gene-id column entirely,
  # so preferring it selected the one file *without* symbols and the first
  # count column silently became the feature ids.
  candidates <- unique(c(raw_files, cuffdiff_files, xlsx_files, norm_files, any_files))
  # Exclude unsupported formats
  candidates <- candidates[!grepl("\\.(mtx|rds|fa|fasta|robj)(\\.(gz|bz2))?$",
                                  candidates, ignore.case = TRUE)]
  candidates <- candidates[!grepl("circRNA|circMIR|novel_circ|_metadata\\.|_meta\\.",
                                  basename(candidates), ignore.case = TRUE)]
  if (length(candidates) == 0) return(NULL)
  candidates[[1]]
}


# ---------------------------------------------------------------------------
# Supplementary file loaders
# ---------------------------------------------------------------------------

# Pain vocabulary for column-name group assignment (R vectors)
PAIN_TOKENS_R <- c(
  "sni", "cci", "snl", "cfa", "formalin", "capsaicin", "carrageenan",
  "psnl", "spnl", "spnt", "axotomy", "neuropath", "nociceptive",
  "allodynia", "hyperalgesia", "fibromyalgia", "osteoarthritis", "rheumatoid",
  "diabetic", "paclitaxel", "oxaliplatin", "mia", "bone.cancer",
  "injured", "ligated", "lesion", "surgery", "ipsilateral"
)
CTRL_TOKENS_R <- c(
  "sham", "naive", "naïve", "control", "ctrl", "ctr", "ctl", "healthy",
  "untreated", "vehicle", "saline", "pbs", "baseline", "normal",
  "uninjured", "mock", "wild.type", "wildtype"
)

# Single-letter abbreviations used in specific count matrices (e.g. GSE186237 S-1 = sham)
CTRL_PREFIX_R <- c("^s[-_.]", "^con[-_.]", "^veh[-_.]")
CASE_PREFIX_R <- c("^t[-_.]", "^trt[-_.]")  # treatment prefix

#' Classify sample column names using pain/control vocabulary.
#' Returns a named character vector: "case" | "control" | NA per column.
classify_cols_by_vocab <- function(col_names) {
  result <- rep(NA_character_, length(col_names))
  lo <- tolower(col_names)

  # Allow digits after tokens: lookbehind/lookahead only blocks letters
  for (tok in PAIN_TOKENS_R) {
    pat <- paste0("(?<![a-z])", tok, "(?![a-z])")
    hits <- grepl(pat, lo, perl = TRUE) & is.na(result)
    result[hits] <- "case"
  }
  for (tok in CTRL_TOKENS_R) {
    pat <- paste0("(?<![a-z])", tok, "(?![a-z])")
    hits <- grepl(pat, lo, perl = TRUE) & is.na(result)
    result[hits] <- "control"
  }

  # Single-character prefix patterns (applied last, only to unassigned columns)
  for (pat in CTRL_PREFIX_R) {
    hits <- grepl(pat, lo, perl = TRUE) & is.na(result)
    result[hits] <- "control"
  }
  for (pat in CASE_PREFIX_R) {
    hits <- grepl(pat, lo, perl = TRUE) & is.na(result)
    result[hits] <- "case"
  }

  setNames(result, col_names)
}


#' Unified supplementary count/expression matrix loader.
#'
#' Handles .txt.gz, .tsv.gz, .csv.gz (tab- or comma-separated) and .xlsx/.xls.
#' Infers sample groups from column names via pain vocabulary; falls back to
#' pheno_used title field if column-name classification fails.
#'
#' Chooses DESeq2 for integer counts, limma-trend for continuous (TPM/FPKM).
#' Assign groups by identity between matrix columns and phenotype rownames.
#'
#' 05d writes {GSE}_raw_count_merged.txt.gz with GSM accessions as column names
#' for some studies and descriptive labels (naive_control_1) for others. GSM
#' columns are the same identifier space as rownames(pData(eset)), so the
#' curated per-sample annotation joins to them directly -- but
#' classify_cols_by_vocab() cannot read a GSM and the title fallback compares
#' free text, so four studies reported "Cannot assign groups (vocab: 0/0;
#' title match: 0/0)" while holding a perfectly good GSM-keyed annotation.
#'
#' Declines (all NA) when the columns are not in the phenotype's space, which
#' is what lets the vocabulary path handle the descriptive-label studies.
#'
#' @param col_names  Column names of the expression matrix.
#' @param pheno_used Phenotype rows with a non-NA `group`, rownames = sample ids.
#' @return Character vector parallel to col_names, named by it; NA where unjoined.
groups_from_sample_ids <- function(col_names, pheno_used) {
  out <- rep(NA_character_, length(col_names))
  names(out) <- col_names
  if (is.null(pheno_used) || nrow(pheno_used) == 0) return(out)
  idx <- match(col_names, rownames(pheno_used))
  out[!is.na(idx)] <- as.character(pheno_used$group[idx[!is.na(idx)]])
  out
}


#' Join matrix column names to phenotype rows through per-sample metadata.
#'
#' A submitter's count matrix names its columns after the library or the tube:
#' `EDCV009_Th1-17_S5`, `7780.EE.0002`, `G1_42L`, `HC-179_FPKM`. None of those
#' is a GSM accession, and none carries pain vocabulary, so both earlier routes
#' decline. GEO holds the same string in the series matrix's own
#' `!Sample_description` -- sometimes bare, sometimes labelled ("Library name:
#' G1_42L", "Column name in counts.txt: 7780.EE.0002") -- so the column can be
#' joined to a sample and the group taken from the curated annotation. Nothing
#' here interprets what a label means; that stays with `{acc}_groups.csv`.
#'
#' Exact matches are taken first and containment only for what is left, so a
#' label that is also a prefix of another (`S1` beside `S1_rep2`) resolves to
#' itself instead of being ambiguous. A column matching more than one sample is
#' left NA rather than assigned to one of them.
#'
#' @param col_names Character vector of matrix column names.
#' @param pheno     Phenotype data.frame with a `group` column, GSM rownames,
#'   and any of `description`, `title`, `supplementary_file` (GEOquery numbers
#'   repeated fields `description.1`, `description.2`, ...).
#' @return Named character vector: "case" | "control" | NA per column.
groups_from_sample_metadata <- function(col_names, pheno,
                                        fields = c("description", "title",
                                                   "supplementary_file")) {
  out <- rep(NA_character_, length(col_names))
  names(out) <- col_names
  if (is.null(pheno) || nrow(pheno) == 0 || !"group" %in% colnames(pheno)) return(out)

  cols <- intersect(fields, colnames(pheno))
  cols <- unique(c(cols, grep("^(description|title)[._][0-9]+$", colnames(pheno),
                              value = TRUE)))
  if (length(cols) == 0) return(out)

  normalise <- function(x) {
    x <- as.character(x)
    x <- ifelse(is.na(x), "", x)
    # "Library name: G1_42L" / "Column name in counts.txt: 7780.EE.0002"
    stripped <- sub("^.{0,60}?:[[:space:]]*", "", x)
    both <- unique(c(x, stripped, basename(stripped)))
    both <- sub("(\\.(gz|txt|csv|tsv|xlsx?))+$", "", both, ignore.case = TRUE)
    tolower(trimws(both[nzchar(both)]))
  }

  # One candidate-string set per sample, over every metadata field present.
  candidates <- lapply(seq_len(nrow(pheno)), function(i) {
    unique(unlist(lapply(cols, function(cn) normalise(pheno[i, cn]))))
  })

  keys <- tolower(trimws(sub("(\\.(gz|txt|csv|tsv|xlsx?))+$", "", col_names,
                            ignore.case = TRUE)))

  assign_from <- function(hits) {
    if (length(hits) != 1) return(NA_character_)
    as.character(pheno$group[hits])
  }

  for (j in seq_along(keys)) {
    if (!nzchar(keys[j])) next
    hits <- which(vapply(candidates, function(cd) keys[j] %in% cd, logical(1)))
    out[j] <- assign_from(hits)
  }

  # Containment only for columns still unresolved, and only for strings long
  # enough that a substring test means something.
  for (j in seq_along(keys)) {
    if (!is.na(out[j]) || nchar(keys[j]) < 3) next
    hits <- which(vapply(candidates, function(cd) {
      cd <- cd[nchar(cd) >= 3]
      length(cd) > 0 &&
        any(vapply(cd, function(s) grepl(s, keys[j], fixed = TRUE) ||
                                   grepl(keys[j], s, fixed = TRUE), logical(1)))
    }, logical(1)))
    out[j] <- assign_from(hits)
  }

  out
}


#' @param curated TRUE when {acc}_groups.csv declared exclusions. The
#'   vocabulary route re-derives labels from column names and so cannot see an
#'   exclusion: GSE217932's curated 4 v 4 became 8 v 4 again because the
#'   THP-treated arm is named like a case. The title and metadata routes carry
#'   `pheno_used$group`, which is the curated label, so when a curator has
#'   excluded samples the vocabulary route is skipped in their favour.
run_deseq2_from_matrix <- function(file_path, pheno_used, accession, cfg,
                                   curated = FALSE) {
  ext <- tolower(file_path)

  # Detect Cuffdiff gene_expression.tsv output and short-circuit
  if (grepl("gene_expression\\.tsv(\\.gz)?$|cuffdiff|diff_exp", basename(file_path),
            ignore.case = TRUE)) {
    res <- run_cuffdiff_output(file_path, accession, cfg)
    if (!is.null(res)) return(res)
  }

  # ---- Read matrix ----
  mat <- tryCatch({
    if (grepl("\\.xlsx?$", ext)) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        message("  readxl not available"); return(NULL)
      }
      suppressMessages(as.data.frame(readxl::read_excel(file_path)))
    } else if (grepl("\\.csv(\\.gz)?$", ext)) {
      read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(e) { message("  Read error: ", e$message); NULL })

  # A single column means the separator was wrong, not that the study has one
  # sample: GSE135080's counts are space-aligned, so the tab reader returned
  # the whole line as one field and the study failed with "Row/ID length
  # mismatch: 0 vs 3858". sep = "" splits on any run of whitespace.
  if (!is.null(mat) && ncol(mat) < 2) {
    alt <- tryCatch(
      read.table(file_path, header = TRUE, sep = "", check.names = FALSE,
                 stringsAsFactors = FALSE, comment.char = ""),
      error = function(e) NULL)
    if (!is.null(alt) && ncol(alt) >= 2) {
      message("  Re-read as whitespace-separated: ", ncol(alt), " columns")
      mat <- alt
    }
  }
  if (is.null(mat) || nrow(mat) == 0) return(NULL)

  # ---- Extract gene IDs and count columns ----
  # If R auto-detected row.names from a leading-comma format, use them directly.
  has_auto_rownames <- !all(rownames(mat) == as.character(seq_len(nrow(mat))))
  if (has_auto_rownames) {
    gene_ids <- make.unique(as.character(rownames(mat)))
  } else {
    gene_ids <- make.unique(as.character(mat[[1]]))
    mat      <- mat[, -1, drop = FALSE]
  }

  # Drop rows with NA/empty gene IDs
  valid_id <- !is.na(gene_ids) & nchar(gene_ids) > 0
  gene_ids <- gene_ids[valid_id]
  mat      <- mat[valid_id, , drop = FALSE]

  # ---- Numeric conversion ----
  # check.names = FALSE, or as.data.frame() runs make.names() over the column
  # names and "con-1" becomes "con.1", "HC-186_Raw.Read.Count" becomes
  # "HC.186_Raw.Read.Count". Every route that joins a column to a sample then
  # compares a mangled string against GEO's unmangled one, so the match fails
  # for exactly the studies whose libraries are named with a hyphen.
  mat_num <- suppressWarnings(
    as.data.frame(lapply(mat, function(x) as.numeric(as.character(x))),
                  check.names = FALSE)
  )
  # Drop annotation-named columns (gene_name, chr, start, end, strand, length, etc.)
  # and any remaining non-numeric columns
  annot_pat <- paste0(
    "^(gene_name|gene_chr|gene_start|gene_end|gene_strand|gene_length|gene_type|gene_biotype",
    "|chr|seqname|chromosome|start|end|strand|length|biotype",
    "|feature|source|frame|attribute|description|entrez|symbol",
    "|ensembl_gene|external|gene_size|width|gc_content|transcript)"
  )
  annot_cols <- grepl(annot_pat, colnames(mat_num), ignore.case = TRUE)
  numeric_cols <- vapply(mat_num, function(x) sum(!is.na(x)) > 0, logical(1))
  keep_cols <- numeric_cols & !annot_cols
  mat_num <- mat_num[, keep_cols, drop = FALSE]
  if (nrow(mat_num) != length(gene_ids)) {
    message("  Row/ID length mismatch: ", nrow(mat_num), " vs ", length(gene_ids), " — skipping")
    return(NULL)
  }
  rownames(mat_num) <- gene_ids
  # Drop all-NA rows
  mat_num <- mat_num[rowSums(!is.na(mat_num)) > 0, , drop = FALSE]
  if (nrow(mat_num) == 0) { message("  Matrix empty after coercion"); return(NULL) }

  # ---- Assign groups ----
  # Identity first: it is the curated per-sample annotation, and it declines
  # cleanly (all NA) when the columns are not sample ids, so the vocabulary
  # path still handles matrices whose columns are descriptive labels.
  col_groups <- groups_from_sample_ids(colnames(mat_num), pheno_used)
  if (sum(col_groups == "case", na.rm = TRUE) < cfg$da$min_samples_per_group ||
      sum(col_groups == "control", na.rm = TRUE) < cfg$da$min_samples_per_group) {
    if (isTRUE(curated)) {
      # Leave col_groups all-NA so the branch below falls through to the title
      # and metadata routes, which honour the curated labels.
      message("  Curated annotation excludes samples — skipping the ",
              "vocabulary route, which cannot see an exclusion")
      col_groups <- rep(NA_character_, ncol(mat_num))
      names(col_groups) <- colnames(mat_num)
    } else {
      col_groups <- classify_cols_by_vocab(colnames(mat_num))
    }
  } else {
    message(sprintf("  Sample-id join: %d case, %d control",
                    sum(col_groups == "case", na.rm = TRUE),
                    sum(col_groups == "control", na.rm = TRUE)))
  }
  n_case_col <- sum(col_groups == "case",    na.rm = TRUE)
  n_ctrl_col <- sum(col_groups == "control", na.rm = TRUE)

  if (n_case_col >= cfg$da$min_samples_per_group &&
      n_ctrl_col >= cfg$da$min_samples_per_group) {
    # Use vocabulary-based assignment
    usable_cols <- names(col_groups)[!is.na(col_groups)]
    mat_use     <- mat_num[, usable_cols, drop = FALSE]
    group_vec   <- col_groups[usable_cols]
    n_case  <- n_case_col
    n_ctrl  <- n_ctrl_col
    message(sprintf("  Column vocab: %d case, %d control", n_case, n_ctrl))
  } else {
    # Fall back: match column names to pheno_used sample titles, then to any
    # per-sample metadata field. The title route stays first and unchanged, so
    # no study it already resolves can be reassigned by the wider one.
    n_case <- 0L
    n_ctrl <- 0L
    matched_cols <- rep(NA_character_, ncol(mat_num))
    if ("title" %in% colnames(pheno_used)) {
      titles <- tolower(pheno_used$title)
      cols_lo <- tolower(colnames(mat_num))
      # Containment alone is ambiguous whenever one label is a substring of
      # another, and the commonest case in this corpus is the plainest:
      # "male" is a substring of "female". Every Male column of GSE197233
      # therefore matched two titles -- its own and the Female one -- and was
      # dropped, halving the study silently. So an ambiguous containment match
      # is re-tried at a word boundary, which "female.sni1.acc" fails for the
      # column "male.sni1" and "male.sni1.acc" passes.
      at_boundary <- function(nm, t) {
        starts <- gregexpr(nm, t, fixed = TRUE)[[1]]
        if (starts[1] == -1) return(FALSE)
        any(starts == 1 | !grepl("[a-z0-9]", substr(t, starts - 1, starts - 1)))
      }
      matched_cols <- vapply(cols_lo, function(nm) {
        # Exact first, the idiom groups_from_sample_metadata() already uses.
        # Without it a column named "mn_1" is ambiguous against the titles
        # "mn_1" and "mn_10" -- both contain it, both at a word boundary -- so
        # it is dropped, and a study whose samples are numbered past nine loses
        # the single-digit ones. GSE250152 lost 18 of 33 that way and
        # GSE186505 two of twenty.
        exact <- which(titles == nm)
        if (length(exact) == 1) return(pheno_used$group[exact])
        idx <- which(sapply(titles, function(t) grepl(t, nm, fixed = TRUE) |
                              grepl(nm, t, fixed = TRUE)))
        if (length(idx) > 1) {
          idx <- idx[sapply(titles[idx], function(t) at_boundary(nm, t))]
        }
        if (length(idx) == 1) pheno_used$group[idx] else NA_character_
      }, character(1))
      n_case <- sum(matched_cols == "case",    na.rm = TRUE)
      n_ctrl <- sum(matched_cols == "control", na.rm = TRUE)
    }

    if (n_case >= cfg$da$min_samples_per_group &&
        n_ctrl >= cfg$da$min_samples_per_group) {
      usable_cols <- colnames(mat_num)[!is.na(matched_cols)]
      mat_use     <- mat_num[, usable_cols, drop = FALSE]
      group_vec   <- matched_cols[!is.na(matched_cols)]
      message(sprintf("  Title match: %d case, %d control", n_case, n_ctrl))
    } else {
      meta_cols <- groups_from_sample_metadata(colnames(mat_num), pheno_used)
      n_case_meta <- sum(meta_cols == "case",    na.rm = TRUE)
      n_ctrl_meta <- sum(meta_cols == "control", na.rm = TRUE)
      if (n_case_meta >= cfg$da$min_samples_per_group &&
          n_ctrl_meta >= cfg$da$min_samples_per_group) {
        usable_cols <- colnames(mat_num)[!is.na(meta_cols)]
        mat_use     <- mat_num[, usable_cols, drop = FALSE]
        group_vec   <- meta_cols[!is.na(meta_cols)]
        n_case <- n_case_meta
        n_ctrl <- n_ctrl_meta
        message(sprintf("  Metadata join: %d case, %d control", n_case, n_ctrl))
      } else {
        message(sprintf(paste0("  Cannot assign groups (vocab: %d/%d; title match: ",
                               "%d/%d; metadata join: %d/%d)"),
                        n_case_col, n_ctrl_col, n_case, n_ctrl,
                        n_case_meta, n_ctrl_meta))
        return(NULL)
      }
    }
  }

  # ---- Decide: integer counts (DESeq2) vs continuous (limma-trend) ----
  vals <- as.numeric(as.matrix(mat_use))
  is_integer_counts <- all(vals == floor(vals), na.rm = TRUE) && max(vals, na.rm = TRUE) > 50

  if (is_integer_counts) {
    int_mat <- matrix(as.integer(vals), nrow = nrow(mat_use), ncol = ncol(mat_use),
                      dimnames = dimnames(mat_use))
    int_mat <- int_mat[rowSums(int_mat, na.rm = TRUE) >= cfg$da$min_count, , drop = FALSE]
    if (nrow(int_mat) == 0) { message("  No features pass min_count"); return(NULL) }

    col_data <- data.frame(group = factor(group_vec, levels = c("control", "case")),
                           row.names = colnames(int_mat))
    dds <- tryCatch(
      DESeqDataSetFromMatrix(countData = int_mat, colData = col_data, design = ~group),
      error = function(e) { message("  DESeq2 setup: ", e$message); NULL }
    )
    if (is.null(dds)) return(NULL)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("group", "case", "control"),
                   alpha = cfg$da$padj_threshold)
    # Use unshrunken MLE for meta-analysis: apeglm posterior SD is too small
    # for null genes and would dominate inverse-variance weighting across studies
    res_df <- as.data.frame(res)
    res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$lfcSE) &
                     res_df$lfcSE > 0, , drop = FALSE]
    message(sprintf("  DESeq2 suppl: %d features", nrow(res_df)))
    build_effect_df(rownames(res_df), res_df$log2FoldChange, res_df$lfcSE,
                    res_df$pvalue, res_df$padj,
                    n_case, n_ctrl, accession, "transcriptomics",
                    detect_id_space_r(rownames(int_mat)), "log2FC")
  } else {
    # TPM / FPKM / normalized — use limma-trend on log2
    log_mat <- log2(mat_use + 1)
    group_fac <- factor(group_vec, levels = c("control", "case"))
    design <- model.matrix(~group_fac)
    fit  <- limma::lmFit(log_mat, design)
    fit2 <- limma::eBayes(fit, trend = TRUE)
    tt   <- limma::topTable(fit2, coef = 2, number = Inf, sort.by = "none")
    se_approx <- se_from_t(tt$logFC, tt$t)
    message(sprintf("  limma-trend suppl: %d features", nrow(tt)))
    build_effect_df(rownames(tt), tt$logFC, se_approx,
                    tt$P.Value, tt$adj.P.Val,
                    n_case, n_ctrl, accession, "transcriptomics",
                    detect_id_space_r(rownames(log_mat)), "log2FC")
  }
}


#' Detect gene ID namespace from a vector of IDs.
detect_id_space_r <- function(ids) {
  ids10 <- head(ids[nchar(ids) > 3], 10)
  if (any(grepl("^ENSMUSG", ids10))) return("mouse_ensembl_gene")
  if (any(grepl("^ENSRNOG", ids10))) return("rat_ensembl_gene")
  if (any(grepl("^ENSG", ids10)))    return("human_ensembl_gene")
  if (any(grepl("^ENSRNOT", ids10))) return("rat_ensembl_transcript")
  if (any(grepl("^ENSMUST", ids10))) return("mouse_ensembl_transcript")
  "symbol"
}

#' limma from CEL files (Affymetrix microarray).
run_limma_from_cel <- function(cel_files, pheno_used, n_case, n_ctrl, accession, cfg) {
  if (!requireNamespace("oligo", quietly = TRUE)) {
    message("  oligo package not available — install via BiocManager::install('oligo')")
    return(NULL)
  }
  library(oligo)

  # Keep only CEL files for usable samples (match by sample title embedded in filename)
  # cel_files are like "GSM4420406_cont1_ACC_MTA-1_0_.CEL.gz"
  usable_gsm <- rownames(pheno_used)
  gsm_in_name <- sapply(cel_files, function(f) {
    m <- regmatches(f, regexpr("GSM[0-9]+", f))
    if (length(m)) m else NA
  })
  keep_cel <- cel_files[!is.na(gsm_in_name) & gsm_in_name %in% usable_gsm]

  if (length(keep_cel) < 2) {
    message("  Too few CEL files matched to usable samples")
    return(NULL)
  }

  raw <- tryCatch(read.celfiles(keep_cel), error = function(e) {
    message("  CEL read error: ", e$message); NULL
  })
  if (is.null(raw)) return(NULL)

  eset_norm <- rma(raw, background = TRUE, normalize = TRUE, target = "core")
  expr <- exprs(eset_norm)

  # Align sample order
  gsm_order <- gsm_in_name[match(colnames(expr), basename(keep_cel))]
  colnames(expr) <- gsm_order
  pheno_used <- pheno_used[colnames(expr), , drop = FALSE]

  group  <- factor(pheno_used$group, levels = c("control", "case"))
  design <- model.matrix(~group)
  fit    <- lmFit(expr, design)
  fit2   <- eBayes(fit)
  tt     <- topTable(fit2, coef = "groupcase", number = Inf, sort.by = "none")

  df <- data.frame(
    feature_id = rownames(tt),
    effect_size = tt$logFC,
    se          = se_from_t(tt$logFC, tt$t),
    pval        = tt$P.Value,
    padj        = tt$adj.P.Val,
    n_case      = n_case, n_control = n_ctrl,
    study_id    = accession, modality = "transcriptomics",
    id_space    = "HGNC", effect_unit = "log2FC",
    stringsAsFactors = FALSE
  )
  message(sprintf("  limma CEL complete: %d features", nrow(df)))
  df
}


# ---------------------------------------------------------------------------
# DESeq2 (RNA-seq)
# ---------------------------------------------------------------------------

run_deseq2 <- function(eset, pheno_used, n_case, n_ctrl, accession, cfg) {
  counts <- exprs(eset)[, rownames(pheno_used)]

  # Keep integer counts only (expression matrix may be TPM/FPKM — flag and skip)
  if (any(counts != floor(counts), na.rm = TRUE)) {
    message("  Non-integer expression values — this study may use TPM/FPKM. Using limma-trend.")
    return(run_limma_trend(eset, pheno_used, n_case, n_ctrl, accession, cfg))
  }

  counts <- counts[rowSums(counts) >= cfg$da$min_count, , drop = FALSE]
  col_data <- data.frame(group = factor(pheno_used$group, levels = c("control", "case")),
                         row.names = rownames(pheno_used))

  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~group),
    error = function(e) { message("  DESeq2 setup failed: ", e$message); NULL }
  )
  if (is.null(dds)) return(NULL)

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "case", "control"),
                 alpha = cfg$da$padj_threshold)

  # Use unshrunken MLE for meta-analysis (apeglm posterior SD too small for null genes)
  df <- as.data.frame(res) %>%
    filter(!is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0) %>%
    transmute(feature_id = rownames(.), effect_size = log2FoldChange, se = lfcSE,
              pval = pvalue, padj = padj)

  build_effect_df(df$feature_id, df$effect_size, df$se, df$pval, df$padj,
                  n_case, n_ctrl, accession, "transcriptomics", "HGNC", "log2FC")
}


# ---------------------------------------------------------------------------
# limma (microarray)
# ---------------------------------------------------------------------------

run_limma_array <- function(eset, pheno_used, n_case, n_ctrl, accession, cfg) {
  expr <- exprs(eset)[, rownames(pheno_used)]

  # GEO series matrices sometimes store linear-scale intensities. See se_utils.R.
  expr <- maybe_log2(expr)

  group <- factor(pheno_used$group, levels = c("control", "case"))
  design <- model.matrix(~group)

  fit  <- lmFit(expr, design)
  fit2 <- eBayes(fit)
  tt   <- topTable(fit2, coef = "groupcase", number = Inf, sort.by = "none")

  build_effect_df(rownames(tt), tt$logFC, se_from_t(tt$logFC, tt$t),
                  tt$P.Value, tt$adj.P.Val,
                  n_case, n_ctrl, accession, "transcriptomics", "HGNC", "log2FC")
}


# ---------------------------------------------------------------------------
# limma-trend (for TPM/FPKM data from RNA-seq without raw counts)
# ---------------------------------------------------------------------------

run_limma_trend <- function(eset, pheno_used, n_case, n_ctrl, accession, cfg) {
  expr <- log2(exprs(eset)[, rownames(pheno_used)] + 1)
  group <- factor(pheno_used$group, levels = c("control", "case"))
  design <- model.matrix(~group)

  fit  <- lmFit(expr, design)
  fit2 <- eBayes(fit, trend = TRUE)
  tt   <- topTable(fit2, coef = "groupcase", number = Inf, sort.by = "none")

  # SE approximated from logFC / t-statistic
  se_approx <- se_from_t(tt$logFC, tt$t)

  build_effect_df(rownames(tt), tt$logFC, se_approx,
                  tt$P.Value, tt$adj.P.Val,
                  n_case, n_ctrl, accession, "transcriptomics", "HGNC", "log2FC")
}
