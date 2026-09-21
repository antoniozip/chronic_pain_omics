#!/usr/bin/env Rscript
#' Normalize per-study feature IDs to gene symbols.
#'
#' Reads effects CSVs from results/per_study/{modality}/ and maps
#' platform-specific IDs (Ensembl transcript IDs, Affymetrix probeset IDs)
#' to gene symbols, then aggregates to gene level (min p-value per gene).
#'
#' Supported mappings:
#'   ENSRNOT* → rat gene symbol (biomaRt rnorvegicus_gene_ensembl)
#'   TC*.mm.* → mouse gene symbol (mta10transcriptcluster.db)
#'   TC*.hg.* → human gene symbol (skip, already has HGNC space)
#'   ENSMUST* → mouse gene symbol (biomaRt mmusculus_gene_ensembl)
#'
#' Output: results/per_study/{modality}/{accession}_effects_normalized.csv
#'   and a combined results/per_study/{modality}/effects_normalized_all.csv
#'
#' Usage:
#'   Rscript pipeline/05b_normalize_ids.R --modality transcriptomics

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
REPO_ROOT <- normalizePath(file.path(dirname(script_path), ".."))

source(file.path(REPO_ROOT, "R", "meta_utils.R"))
source(file.path(REPO_ROOT, "R", "symbol_utils.R"))

option_list <- list(
  make_option("--modality",     type = "character", default = "transcriptomics"),
  make_option("--per-study-dir", type = "character",
              default = file.path(REPO_ROOT, "results", "per_study")),
  make_option("--interim-dir",  type = "character",
              default = file.path(REPO_ROOT, "data", "interim")),
  make_option("--cache-dir",    type = "character",
              default = file.path(REPO_ROOT, "data", "raw", "geo_cache"))
)
opt <- parse_args(OptionParser(option_list = option_list))

per_study_dir <- opt[["per-study-dir"]]
interim_dir   <- opt[["interim-dir"]]
cache_dir     <- opt[["cache-dir"]]
modality      <- opt$modality

message(sprintf("\n========== NORMALIZE IDs: %s ==========", toupper(modality)))


#' Species per accession, read from the harmonized manifest.
#'
#' Entrez ids name no organism, so mapping them needs the study's species from
#' somewhere; this step previously knew only the GPL. Derived study ids
#' (GSE241361_DRG) are absent from the manifest, so fall back to the parent
#' accession by stripping the suffix, as 06g does.
load_species_map <- function(interim_dir, modality) {
  path <- file.path(interim_dir, modality, "harmonized.jsonl")
  if (!file.exists(path)) {
    message("No harmonized manifest at ", path, " — species-dependent maps unavailable.")
    return(list())
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  out <- list()
  for (ln in lines) {
    rec <- jsonlite::fromJSON(ln, simplifyVector = TRUE)
    out[[rec$accession]] <- as.character(rec$species_canonical)
  }
  out
}

SPECIES_MAP <- load_species_map(interim_dir, modality)

species_for <- function(accession) {
  if (!is.null(SPECIES_MAP[[accession]])) return(SPECIES_MAP[[accession]])
  parent <- sub("_[A-Za-z_]+$", "", accession)
  if (!is.null(SPECIES_MAP[[parent]])) return(SPECIES_MAP[[parent]])
  NA_character_
}

# ---------------------------------------------------------------------------
# Helper: detect ID space from feature IDs
# ---------------------------------------------------------------------------

#' Strip the NCBI GFF `gene-` prefix that some pipelines carry into a count
#' matrix (GSE236592: `gene-Dgkh`).
#'
#' Only `gene-`. `rna-` and `cds-` name a transcript and a CDS, so stripping
#' them yields an accession rather than a symbol and they must fall through to
#' the normal detection. A proportion, not `any`, for the reason the rest of
#' this file uses proportions: a handful of prefixed rows in an otherwise
#' symbol-keyed table is not a prefixed table.
#'
#' @param ids Character vector of feature ids.
#' @param min_share Minimum fraction carrying the prefix before it is stripped.
#' @return The ids, with `gene-` removed if the table is predominantly prefixed.
strip_gene_prefix <- function(ids, min_share = 0.95) {
  present <- ids[!is.na(ids) & nzchar(ids)]
  if (length(present) == 0) return(ids)
  if (mean(grepl("^gene-", present)) < min_share) return(ids)
  sub("^gene-", "", ids)
}


#' Share of a study's ids that must lie in a space for it to be routed there.
#' A strict majority, so a stray minority cannot route the whole study: that is
#' the mistake in both directions, and both directions are silent. Strict
#' because a tie is not a majority -- half digits and half symbols is a mixed
#' table, not an Entrez study.
#'
#' Below this, the study is `unknown` and keeps its native ids. Above it but
#' below `ID_SPACE_UNANIMOUS`, the routing happens and says how many ids it will
#' drop, because a non-unanimous majority is a mixed table and no mapper here
#' handles one.
ID_SPACE_MIN_SHARE  <- 0.5
ID_SPACE_UNANIMOUS  <- 0.95

#' Share of `ids` matching `pattern`, over the non-empty ids.
id_share <- function(ids, pattern) {
  present <- ids[!is.na(ids) & nzchar(ids)]
  if (length(present) == 0) return(0)
  mean(grepl(pattern, present))
}

#: Ordered space -> pattern. The five Ensembl prefixes are mutually exclusive,
#: so at most one can hold a majority; order still matters for the digit and
#: accession patterns further down.
ID_SPACE_PATTERNS <- c(
  rat_ensembl_transcript   = "^ENSRNOT",
  mouse_mta10              = "^TC.*\\.mm\\.",
  mouse_ensembl_transcript = "^ENSMUST",
  mouse_ensembl_gene       = "^ENSMUSG",
  rat_ensembl_gene         = "^ENSRNOG",
  # ENSG is safe against the rodent prefixes above, which are ENSM.. / ENSR..
  human_ensembl_gene       = "^ENSG[0-9]",
  # Human Ensembl *transcripts*, reported separately from genes because
  # org.Hs.eg.db carries only ~39k ENSEMBLTRANS keys and maps almost none of a
  # full transcript-level table. Naming the space is the point: "unknown" is
  # indistinguishable in the output from a study whose ids were already
  # symbols, and GSE153739/GSE153740 sat in the pool as 299,425 unjoinable
  # transcript ids for exactly that reason.
  human_ensembl_transcript = "^ENST[0-9]",
  # Bare Entrez ids. Array-address ids are also digits, so a study mixing
  # digits with symbols is not an Entrez study; map_entrez() additionally
  # refuses a poor match rate, because an array address can collide with a real
  # Entrez id and map to an unrelated gene.
  entrez                   = "^[0-9]+$",
  # RefSeq transcripts, versioned or not (GSE232713: NM_000015.2). The version
  # separator is either form: GSE193928 deposits NM_022892_2, and requiring a
  # dot put it at 88%, so a study that is 100% RefSeq read as unknown.
  refseq                   = "^[NX][MR]_[0-9]+([._][0-9]+)?$"
)

#' Which identifier space a study's feature ids are in.
#'
#' Every test is a share of the study's ids, never `any()` or `all()`. Both
#' extremes have already shipped defects here and neither raises:
#'
#'   * `any()` lets a stray minority route the whole study. GSE221921 is 99.8%
#'     HGNC symbols with 42 ENSG rows out of 21,914; GSE198608 is 0.72%
#'     ENSRNOG, and being sent down the rat Ensembl route took it from 21,170
#'     features to **2**. GSE289659, at 13.7%, went from 18,142 to 378.
#'   * `all()` lets a stray minority veto the mapping. GSE117526 carries 6
#'     RefSeq accessions among 14,346 Entrez ids, so 0.04% of the study
#'     blocked it.
#'
#' The threshold is a majority rather than 0.95 because two studies in the
#' corpus are legitimately 90.5% (GSE147216, MTA-1.0) and 94.3% (GSE245768,
#' rat Ensembl) that space, and no study anywhere in the corpus sits between
#' 20% and 90% -- so a majority separates the real cases from the strays with
#' the widest margin available.
detect_id_space <- function(ids) {
  for (space in names(ID_SPACE_PATTERNS)) {
    share <- id_share(ids, ID_SPACE_PATTERNS[[space]])
    if (share <= ID_SPACE_MIN_SHARE) next
    if (share < ID_SPACE_UNANIMOUS) {
      message(sprintf(paste0("  %s holds %.1f%% of the ids, so it is a mixed ",
                             "table: the %.1f%% that are not will not map"),
                      space, 100 * share, 100 * (1 - share)))
    }
    return(space)
  }
  "unknown"
}


#' Ensembl mart dataset and cache tag for a species, or NULL if unusable.
#'
#' Explicit, and deliberately not the old `if (mouse) mouse else rat`: that
#' two-way form sent every unrecognised species -- human included -- to the rat
#' mart, which returns almost nothing and is indistinguishable in the output
#' from a study with poor annotation. Two of the five studies that reached the
#' meta-analysis in Ensembl space were human.
#'
#' A study recording more than one species (GSE236754) returns NULL rather than
#' silently taking the first: mapping half its features to the wrong organism
#' is worse than declining.
#'
#' @param species Character vector or scalar, e.g. "Homo sapiens".
#' @return Named character vector with `tag` and `dataset`, or NULL.
ensembl_dataset_for <- function(species) {
  s <- tolower(paste(species[!is.na(species)], collapse = " "))
  if (!nzchar(s)) return(NULL)
  hits <- c(
    mouse = grepl("mus|mouse", s),
    rat   = grepl("rattus|\\brat\\b", s),
    human = grepl("sapiens|human", s)
  )
  if (sum(hits) != 1) return(NULL)
  key <- names(hits)[hits]
  switch(key,
    mouse = c(tag = "mouse", dataset = "mmusculus_gene_ensembl"),
    rat   = c(tag = "rat",   dataset = "rnorvegicus_gene_ensembl"),
    human = c(tag = "human", dataset = "hsapiens_gene_ensembl")
  )
}


#' Map bare Entrez gene ids to symbols via the org.*.eg.db for the species.
#'
#' Declines below `min_match` rather than returning a sparse map: bare digits
#' are also how some arrays name their addresses, and an address that happens
#' to be a valid Entrez id maps to an unrelated gene. A low match rate is the
#' signal that the ids are not really Entrez.
#'
#' @param ids       Character vector of feature ids.
#' @param species   Species string(s) for the study.
#' @param min_match Minimum fraction of ids that must resolve to a symbol.
#' @return Named character vector keyed by id, or NULL.
map_entrez <- function(ids, species, min_match = 0.5) {
  d <- ensembl_dataset_for(species)
  if (is.null(d)) {
    message("  Entrez-like ids but species is unusable (", 
            paste(species, collapse = "; "), ") — cannot map")
    return(NULL)
  }
  pkg <- c(mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db",
           human = "org.Hs.eg.db")[[d[["tag"]]]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  ", pkg, " not installed — cannot map Entrez ids")
    return(NULL)
  }
  uniq <- unique(ids)
  sym <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(
      getExportedValue(pkg, pkg), keys = uniq, column = "SYMBOL",
      keytype = "ENTREZID", multiVals = "first")),
    error = function(e) { message("  ", pkg, " lookup failed: ", e$message); NULL })
  if (is.null(sym)) return(NULL)
  matched <- sum(!is.na(sym)) / length(uniq)
  message(sprintf("  Entrez → symbol via %s: %d / %d (%.0f%%)",
                  pkg, sum(!is.na(sym)), length(uniq), 100 * matched))
  if (matched < min_match) {
    message("  Match rate below ", min_match,
            " — treating these as non-Entrez ids and declining")
    return(NULL)
  }
  setNames(as.character(sym), uniq)
}

#' Map RefSeq transcript accessions to symbols via the org.*.eg.db.
#'
#' Versions are stripped before lookup, because org.*.eg.db keys REFSEQ
#' unversioned while submitters deposit either form. Declines below `min_match`
#' for the same reason `map_entrez()` does: a sparse map is the signal that the
#' ids are not what they looked like.
#'
#' @param ids       Character vector of feature ids.
#' @param species   Species string(s) for the study.
#' @param min_match Minimum fraction of ids that must resolve to a symbol.
#' @return Named character vector keyed by the original (possibly versioned)
#'   id, or NULL.
map_refseq <- function(ids, species, min_match = 0.5) {
  d <- ensembl_dataset_for(species)
  if (is.null(d)) {
    message("  RefSeq-like ids but species is unusable (",
            paste(species, collapse = "; "), ") — cannot map")
    return(NULL)
  }
  pkg <- c(mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db",
           human = "org.Hs.eg.db")[[d[["tag"]]]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  ", pkg, " not installed — cannot map RefSeq ids")
    return(NULL)
  }
  uniq  <- unique(ids)
  bare  <- sub("\\.[0-9]+$", "", uniq)
  sym <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(
      getExportedValue(pkg, pkg), keys = unique(bare), column = "SYMBOL",
      keytype = "REFSEQ", multiVals = "first")),
    error = function(e) { message("  ", pkg, " lookup failed: ", e$message); NULL })
  if (is.null(sym)) return(NULL)
  out <- setNames(as.character(sym[bare]), uniq)
  matched <- sum(!is.na(out)) / length(uniq)
  message(sprintf("  RefSeq → symbol via %s: %d / %d (%.0f%%)",
                  pkg, sum(!is.na(out)), length(uniq), 100 * matched))
  if (matched < min_match) {
    message("  Match rate below ", min_match,
            " — treating these as non-RefSeq ids and declining")
    return(NULL)
  }
  out
}


# ---------------------------------------------------------------------------
# Mapping: rat Ensembl transcript → rat gene symbol
# ---------------------------------------------------------------------------

map_rat_ensembl_transcript <- function(ids, interim_dir) {
  cache_path <- file.path(interim_dir, "rat_ensrnot2symbol.csv")
  if (file.exists(cache_path)) {
    message("  Loading cached rat transcript → symbol map")
    cache <- read.csv(cache_path, stringsAsFactors = FALSE)
    base_map <- setNames(cache$external_gene_name, cache$ensembl_transcript_id)
    # Reindex to input IDs (versioned) by stripping version suffix on lookup
    ids_noversion <- sub("\\.[0-9]+$", "", unique(ids))
    return(setNames(base_map[ids_noversion], unique(ids)))
  }

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("  biomaRt not available — cannot map rat transcript IDs")
    return(NULL)
  }

  # Ensembl stores IDs without version suffix; strip .version before querying
  ids_noversion <- sub("\\.[0-9]+$", "", unique(ids))

  message("  Querying biomaRt for rat Ensembl transcript → gene symbol ...")
  # Allow adequate time for Ensembl server
  old_timeout <- getOption("timeout")
  options(timeout = 300)
  on.exit(options(timeout = old_timeout))

  mart <- tryCatch(
    biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl"),
    error = function(e) {
      message("  biomaRt connection failed: ", e$message); NULL
    }
  )
  if (is.null(mart)) return(NULL)

  # Query in batches to avoid timeouts on large ID lists
  batch_size <- 1000
  results <- list()
  for (i in seq(1, length(ids_noversion), batch_size)) {
    batch <- ids_noversion[i:min(i + batch_size - 1, length(ids_noversion))]
    res <- tryCatch(
      biomaRt::getBM(
        attributes = c("ensembl_transcript_id", "external_gene_name"),
        filters    = "ensembl_transcript_id",
        values     = batch,
        mart       = mart
      ),
      error = function(e) { message("  biomaRt batch failed: ", e$message); NULL }
    )
    if (!is.null(res) && nrow(res) > 0) results[[length(results) + 1]] <- res
  }
  if (length(results) == 0) return(NULL)
  res <- do.call(rbind, results)
  res <- res[res$external_gene_name != "", ]

  dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res, cache_path, row.names = FALSE)
  message(sprintf("  Mapped %d / %d transcripts to rat gene symbols", nrow(res), length(ids_noversion)))

  # Build map from original (versioned) IDs to gene symbol
  id_noversion <- sub("\\.[0-9]+$", "", unique(ids))
  symbol_map   <- setNames(res$external_gene_name, res$ensembl_transcript_id)
  setNames(symbol_map[id_noversion], unique(ids))
}

# ---------------------------------------------------------------------------
# Mapping: mouse MTA-1.0 probeset → mouse gene symbol
# ---------------------------------------------------------------------------

map_mouse_mta10 <- function(ids, interim_dir) {
  map_path <- file.path(interim_dir, "mta10_probe2symbol.csv")
  if (!file.exists(map_path)) {
    if (!requireNamespace("mta10transcriptcluster.db", quietly = TRUE)) {
      message("  mta10transcriptcluster.db not available — skipping MTA-1.0 mapping")
      return(NULL)
    }
    message("  Building MTA-1.0 probeset → mouse gene symbol cache ...")
    suppressPackageStartupMessages(library(mta10transcriptcluster.db))
    all_probes <- AnnotationDbi::keys(mta10transcriptcluster.db, keytype = "PROBEID")
    tc_mm <- all_probes[grep("^TC.*mm", all_probes)]
    results <- list()
    batch_size <- 5000
    for (i in seq(1, length(tc_mm), batch_size)) {
      batch <- tc_mm[i:min(i + batch_size - 1, length(tc_mm))]
      res <- tryCatch(
        AnnotationDbi::select(mta10transcriptcluster.db,
          keys    = batch, keytype = "PROBEID",
          columns = c("PROBEID", "SYMBOL")),
        error = function(e) NULL
      )
      if (!is.null(res)) results[[length(results) + 1]] <- res
    }
    full_map <- do.call(rbind, results)
    full_map <- full_map[!is.na(full_map$SYMBOL), ]
    dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(full_map, map_path, row.names = FALSE)
    message(sprintf("  Saved %d probe→symbol mappings", nrow(full_map)))
  } else {
    message("  Loading cached MTA-1.0 → symbol map")
  }

  cache <- read.csv(map_path, stringsAsFactors = FALSE)
  id_map <- setNames(cache$SYMBOL, cache$PROBEID)

  n_mapped <- sum(ids %in% names(id_map))
  message(sprintf("  Mapped %d / %d probesets to mouse gene symbols", n_mapped, length(ids)))
  id_map
}

# ---------------------------------------------------------------------------
# Mapping: Ensembl gene IDs → gene symbols via biomaRt (mouse or rat)
# ---------------------------------------------------------------------------

#' Ensembl gene id -> symbol from an org.*.eg.db, without the network.
#'
#' biomaRt::useMart() answers HTTP 404 on this workstation, so the Ensembl
#' route is unavailable even with connectivity; the org packages carry the same
#' mapping locally. Tried before biomaRt, after the on-disk cache, so the
#' existing mouse and rat caches keep their exact behaviour.
#'
#' Declines below `min_match` so a wrong species or a non-Ensembl id set falls
#' through to biomaRt rather than returning a mostly-empty map.
#'
#' @param ids         Character vector of Ensembl gene ids, versioned or not.
#' @param species_tag One of "mouse", "rat", "human".
#' @param min_match   Minimum fraction that must resolve.
#' @return Named character vector keyed by the original ids, or NULL.
map_ensembl_via_orgdb <- function(ids, species_tag, min_match = 0.3) {
  pkg <- c(mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db",
           human = "org.Hs.eg.db")[[species_tag]]
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  uniq      <- unique(as.character(ids))
  noversion <- sub("\\.[0-9]+$", "", uniq)
  sym <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(
      getExportedValue(pkg, pkg), keys = unique(noversion), column = "SYMBOL",
      keytype = "ENSEMBL", multiVals = "first")),
    error = function(e) NULL)
  if (is.null(sym)) return(NULL)
  out <- setNames(as.character(sym[noversion]), uniq)
  matched <- sum(!is.na(out)) / length(out)
  if (matched < min_match) return(NULL)
  message(sprintf("  Ensembl → symbol via %s: %d / %d (%.0f%%)",
                  pkg, sum(!is.na(out)), length(out), 100 * matched))
  out
}


#' Ensembl transcript accessions to gene symbols.
#'
#' Separate from map_ensembl_gene() because the keytype differs: ENSEMBL is the
#' gene keytype and returns nothing for a transcript accession.
#'
#' Prefers the interim cache, then falls back to the species org package. The
#' fallback rarely succeeds on its own -- org.Hs.eg.db carries ~39k
#' ENSEMBLTRANS keys against >250k human transcripts and resolved 1 of 5
#' sampled GSE153740 ids -- which is why the cache exists. Declining below the
#' match floor is the correct outcome and is reported, because the alternative
#' is writing transcript accessions out as if they were symbols: GSE153739 and
#' GSE153740 reached the meta-analysis that way, contributing 215,320 features
#' that no other study could ever join.
map_ensembl_transcript <- function(ids, species_tag, interim_dir,
                                   min_match = 0.3) {
  # The cache first, exactly as map_rat_ensembl_transcript() does. It is built
  # by scripts/build_transcript_symbol_map.py from MyGene.info, because the two
  # sources this function would otherwise use both fail here: biomaRt answers
  # HTTP 404, and the org packages are thin at transcript level.
  cache_name <- c(mouse = "mouse_ensmust2symbol.csv",
                  rat   = "rat_ensrnot2symbol.csv",
                  human = "human_enst2symbol.csv")[[species_tag]]
  cache_path <- file.path(interim_dir, cache_name)
  if (file.exists(cache_path)) {
    message("  Loading cached ", species_tag, " transcript → symbol map")
    cache <- read.csv(cache_path, stringsAsFactors = FALSE)
    base_map <- setNames(cache$external_gene_name, cache$ensembl_transcript_id)
    uniq <- unique(as.character(ids))
    out  <- setNames(base_map[sub("\\.[0-9]+$", "", uniq)], uniq)
    matched <- sum(!is.na(out)) / length(out)
    message(sprintf("  Transcript → symbol from cache: %d / %d (%.1f%%)",
                    sum(!is.na(out)), length(out), 100 * matched))
    if (matched >= min_match) return(out)
    message("  Below the ", min_match, " match floor — falling through.")
  }

  pkg <- c(mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db",
           human = "org.Hs.eg.db")[[species_tag]]
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  uniq      <- unique(as.character(ids))
  noversion <- sub("\\.[0-9]+$", "", uniq)
  sym <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(
      getExportedValue(pkg, pkg), keys = unique(noversion), column = "SYMBOL",
      keytype = "ENSEMBLTRANS", multiVals = "first")),
    error = function(e) NULL)
  if (is.null(sym)) return(NULL)
  out <- setNames(as.character(sym[noversion]), uniq)
  matched <- sum(!is.na(out)) / length(out)
  message(sprintf("  Ensembl transcript → symbol via %s: %d / %d (%.1f%%)",
                  pkg, sum(!is.na(out)), length(out), 100 * matched))
  if (matched < min_match) {
    message("  Below the ", min_match, " match floor — declining, so the study ",
            "is not written out in transcript space as if it were symbols.")
    return(NULL)
  }
  out
}


map_ensembl_gene <- function(ids, species, interim_dir) {
  d <- ensembl_dataset_for(species)
  if (is.null(d)) {
    message("  Cannot resolve an Ensembl mart for species: ",
            paste(species, collapse = "; "))
    return(NULL)
  }
  species_tag  <- d[["tag"]]
  dataset_name <- d[["dataset"]]
  cache_path   <- file.path(interim_dir, paste0(species_tag, "_ensembl_gene2symbol.csv"))

  if (file.exists(cache_path)) {
    message("  Loading cached ", species_tag, " Ensembl gene → symbol map")
    cache   <- read.csv(cache_path, stringsAsFactors = FALSE)
    id_map  <- setNames(cache$external_gene_name, cache$ensembl_gene_id)
    # Strip version suffix for lookup
    ids_noversion <- sub("\\.[0-9]+$", "", unique(ids))
    return(setNames(id_map[ids_noversion], unique(ids)))
  }

  # Offline first: biomaRt is not reachable here (HTTP 404 from useMart).
  om <- map_ensembl_via_orgdb(ids, species_tag)
  if (!is.null(om)) return(om)

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("  biomaRt not available — cannot map Ensembl gene IDs")
    return(NULL)
  }

  message("  Querying biomaRt (", dataset_name, ") for Ensembl gene → symbol ...")
  old_timeout <- getOption("timeout"); options(timeout = 300)
  on.exit(options(timeout = old_timeout))

  ids_noversion <- sub("\\.[0-9]+$", "", unique(ids))
  mart <- tryCatch(
    biomaRt::useMart("ensembl", dataset = dataset_name),
    error = function(e) { message("  biomaRt failed: ", e$message); NULL }
  )
  if (is.null(mart)) return(NULL)

  batch_size <- 1000
  results <- list()
  for (i in seq(1, length(ids_noversion), batch_size)) {
    batch <- ids_noversion[i:min(i + batch_size - 1, length(ids_noversion))]
    res <- tryCatch(
      biomaRt::getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters    = "ensembl_gene_id",
        values     = batch,
        mart       = mart
      ),
      error = function(e) NULL
    )
    if (!is.null(res) && nrow(res) > 0) results[[length(results) + 1]] <- res
  }
  if (length(results) == 0) return(NULL)
  res <- do.call(rbind, results)
  res <- res[res$external_gene_name != "", ]
  write.csv(res, cache_path, row.names = FALSE)
  message(sprintf("  Mapped %d / %d gene IDs to symbols", nrow(res), length(ids_noversion)))

  id_map <- setNames(res$external_gene_name, res$ensembl_gene_id)
  setNames(id_map[ids_noversion], unique(ids))
}


# ---------------------------------------------------------------------------
# Platform-based dispatch (GEO GPL → Bioconductor annotation .db)
# ---------------------------------------------------------------------------

GPL_TO_ANNO_PKG <- list(
  # Affymetrix Rat arrays
  "GPL85"   = "rgu34a.db",                      # Rat Genome U34A
  "GPL86"   = "rgu34b.db",                      # Rat Genome U34B
  "GPL87"   = "rgu34c.db",                      # Rat Genome U34C
  "GPL341"  = "rae230a.db",                     # Rat Expression 230A
  "GPL1355" = "rat2302.db",                     # Rat Genome 230 2.0
  # Affymetrix Mouse arrays
  "GPL1261" = "mouse4302.db",                   # Mouse Genome 430 2.0
  "GPL6246" = "mogene10sttranscriptcluster.db", # Mouse Gene 1.0 ST
  # GPL6247 was listed here as a "Mouse Gene 1.0 ST variant". It is the
  # Affymetrix *Rat* Gene 1.0 ST array, so the lookup could only ever fail;
  # it now resolves through its GEO platform table's gene_assignment column.
  # GPL6194 and GPL6543 (Affymetrix Rat Exon 1.0 ST) and GPL2050
  # (Rosetta/Merck Rat 25k) were likewise mapped to Illumina *mouse* packages
  # and are removed for the same reason. Checked against
  # !Platform_title in each cached SOFT file on 2026-08-15.
  "GPL16570" = "mogene11sttranscriptcluster.db",# Mouse Gene 1.1 ST
  # Illumina BeadArrays (mouse)
  "GPL6885" = "illuminaMousev2.db",             # MouseRef-8 v2.0
  "GPL6887" = "illuminaMousev2.db"              # MouseWG-6 v2.0 (variant)
)

#' Read the GPL platform accession from a GEO series matrix file.
#' Map feature ids to gene symbols using the cached GEO platform table.
#'
#' The annotation-package route fails for two large classes of array. Some
#' platforms label features by array-address number (`5617701`, `10700001`)
#' while the .db packages key on strings like `10181072_239_rc-S`, so the
#' lookup returns "None of the keys entered are valid keys for 'PROBEID'".
#' Others have no Bioconductor package at all. In both cases GEO's own
#' platform table already carries the mapping and is already cached, because
#' `detect_platform()` downloads it to read one header line.
#'
#' Two table shapes cover the platforms in this corpus:
#'   - an explicit symbol column (`GeneSymbol`, `Gene symbol`, `SYMBOL`, ...)
#'   - Affymetrix `gene_assignment`, which packs several `//`-delimited fields
#'     per row (`NM_012345 // Actb // actin beta // ...`); the symbol is the
#'     second, and rows can carry several assignments separated by `///`.
#'
#' Returns NULL rather than an empty map when nothing resolves, so the caller's
#' existing "keeping native IDs" path still applies.
#' Path to a platform's cached SOFT table, fetching it if it is not there yet.
#'
#' `detect_platform()` reads the GPL accession out of the series matrix and
#' does **not** download anything, so the comment that this table "is already
#' cached" only held for platforms an earlier run happened to fetch. The
#' 2026-08-28 search amendment brought thirteen platforms that no run had
#' touched, and every study on them fell through to native probe ids -- the
#' 2026-08-16 probe-space defect returning through a different door, and just
#' as silent, because keeping native ids is also what a study whose ids are
#' already symbols looks like.
#'
#' @param gpl       Platform accession, e.g. "GPL570".
#' @param cache_dir Directory holding `{gpl}.soft.gz`.
#' @return Path to the cached file, or NA_character_ if it cannot be obtained.
ensure_gpl_soft <- function(gpl, cache_dir) {
  if (is.null(gpl) || is.na(gpl)) return(NA_character_)
  for (p in file.path(cache_dir, paste0(gpl, c(".soft.gz", ".soft", ".annot.gz")))) {
    if (file.exists(p)) return(p)
  }
  if (!requireNamespace("GEOquery", quietly = TRUE)) return(NA_character_)
  message("  Fetching platform table for ", gpl)
  ok <- tryCatch({
    suppressMessages(GEOquery::getGEO(gpl, destdir = cache_dir))
    TRUE
  }, error = function(e) {
    message("  Platform fetch failed for ", gpl, ": ", e$message)
    FALSE
  })
  if (!ok) return(NA_character_)
  for (p in file.path(cache_dir, paste0(gpl, c(".soft.gz", ".soft", ".annot.gz")))) {
    if (file.exists(p)) return(p)
  }
  NA_character_
}


map_via_gpl_table <- function(ids, gpl, cache_dir) {
  if (is.na(gpl) || is.null(gpl)) return(NULL)
  path <- ensure_gpl_soft(gpl, cache_dir)
  if (is.na(path) || !file.exists(path)) return(NULL)

  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  start <- grep("^!platform_table_begin", lines)
  if (!length(start)) return(NULL)
  end <- grep("^!platform_table_end", lines)
  end <- if (length(end)) end[1] - 1L else length(lines)
  body <- lines[(start[1] + 1L):end]
  if (length(body) < 2) return(NULL)

  tbl <- tryCatch(
    utils::read.delim(text = paste(body, collapse = "\n"), quote = "",
                      check.names = FALSE, stringsAsFactors = FALSE,
                      colClasses = "character"),
    error = function(e) NULL)
  if (is.null(tbl) || !nrow(tbl) || !("ID" %in% colnames(tbl))) return(NULL)

  sym_col <- intersect(
    c("GeneSymbol", "Gene symbol", "Gene Symbol", "SYMBOL", "Symbol",
      "gene_symbol", "GENE_SYMBOL",
      # Arraystar platforms (GPL25915) carry the symbol in ORF and have no
      # column from the list above. Last in the list so an explicit symbol
      # column always wins where both exist.
      "ORF"), colnames(tbl))
  if (length(sym_col)) {
    symbols <- trimws(tbl[[sym_col[1]]])
    source_col <- sym_col[1]
  } else if ("gene_assignment" %in% colnames(tbl)) {
    # Take the first assignment's symbol; "---" is Affymetrix for unassigned.
    first <- sub("///.*$", "", tbl[["gene_assignment"]])
    parts <- strsplit(first, "//", fixed = TRUE)
    symbols <- trimws(vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_,
                             character(1)))
    source_col <- "gene_assignment"
  } else {
    return(NULL)
  }

  symbols[symbols %in% c("", "---", "NA", "null")] <- NA_character_
  keep <- !is.na(symbols) & !is.na(tbl$ID) & nzchar(trimws(tbl$ID))
  if (!any(keep)) return(NULL)

  # A named character vector keyed by feature id, matching what every other
  # mapper here returns; aggregate_to_gene() indexes it as `id_map[ids]`, so a
  # data.frame silently fails with "undefined columns selected".
  probe_ids <- trimws(tbl$ID[keep])
  syms <- symbols[keep]
  dup <- duplicated(probe_ids)
  map <- setNames(syms[!dup], probe_ids[!dup])

  hits <- sum(unique(as.character(ids)) %in% names(map))
  if (hits == 0) return(NULL)
  message(sprintf("  Mapped %d / %d ids via %s platform table (%s)",
                  hits, length(unique(as.character(ids))), gpl, source_col))
  map
}


detect_platform <- function(accession, cache_dir) {
  matrix_files <- list.files(
    cache_dir,
    pattern = paste0("^", accession, "(-GPL\\w+)?_series_matrix\\.txt\\.gz$"),
    full.names = TRUE
  )
  if (length(matrix_files) == 0) return(NA_character_)

  con <- gzfile(matrix_files[1], "rt")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_platform_id")) {
      parts <- strsplit(line, "\t")[[1]]
      return(gsub('"', "", parts[2]))
    }
  }
  NA_character_
}


#' Generic probe → symbol mapper using an AnnotationDbi-style .db package.
#' Returns a named character vector keyed by input probe IDs.
map_via_annotation_db <- function(ids, pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message("  Annotation package not available: ", pkg_name)
    return(NULL)
  }
  suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
  anno_db <- get(pkg_name)
  unique_ids <- unique(ids)

  res <- tryCatch(
    AnnotationDbi::select(anno_db, keys = unique_ids,
                          keytype = "PROBEID",
                          columns = c("PROBEID", "SYMBOL")),
    error = function(e) { message("  Lookup failed: ", e$message); NULL }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res <- res[!is.na(res$SYMBOL) & res$SYMBOL != "", ]
  if (nrow(res) == 0) return(NULL)
  # Keep first symbol per probe (probes with multiple gene hits collapse to one)
  res <- res[!duplicated(res$PROBEID), ]
  message(sprintf("  Mapped %d / %d probes via %s",
                  nrow(res), length(unique_ids), pkg_name))
  setNames(res$SYMBOL, res$PROBEID)
}


# ---------------------------------------------------------------------------
# Aggregate to gene level: min pval per gene symbol
# ---------------------------------------------------------------------------

#' Collapse probes to one row per gene, selecting on precision rather than
#' on the result.
#'
#' This used to keep the probe with the smallest p-value, which is selection on
#' the outcome and biases in two directions at once:
#'
#'   - **The p-value stops being a p-value.** The minimum of n draws is
#'     Beta(1, n) under the null, not uniform. GSE53860 carries 226,599 exon
#'     probes over 15,123 genes — 15 per gene — so its "p-values" had an
#'     expectation near 0.06 instead of 0.5, and genes with the most probes won
#'     most often. That is why the SNL stratum jumped from 273 to 1,833
#'     significant features when those studies were first mapped, with
#'     ribosomal genes (Rpl30, Rpl41) at the top: they have the most probes,
#'     not the strongest signal.
#'   - **The effect size is biased away from zero.** `06` pools `effect_size`
#'     and `se`, not `pval`, so picking the most significant probe hands the
#'     meta-analysis the most extreme estimate of the gene. Correcting the
#'     p-value alone would leave this untouched.
#'
#' Selecting the smallest standard error fixes both. Precision depends on the
#' measurement's variance and sample size, not on where the effect landed, so
#' the retained estimate is unbiased for the gene and the p-value keeps its
#' meaning. Ties break on the feature id so the choice is reproducible.
#'
#' Inverse-variance averaging across a gene's probes would use more of the data
#' and is tempting, but probes of one transcript are strongly correlated;
#' treating them as independent would shrink the standard error toward a
#' precision the study does not have, which is the pseudo-replication this
#' project rejects everywhere else.
aggregate_to_gene <- function(df, id_map) {
  df$gene_symbol <- id_map[df$feature_id]
  df_mapped <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", ]
  if (nrow(df_mapped) == 0) return(NULL)

  # Probes with an unusable SE cannot be ranked on precision and would sort to
  # the front of an ascending order; drop them unless a gene has nothing else.
  df_mapped$.se_rank <- ifelse(is.na(df_mapped$se) | df_mapped$se <= 0,
                               Inf, df_mapped$se)

  df_mapped %>%
    group_by(gene_symbol) %>%
    arrange(.se_rank, feature_id, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(-.se_rank) %>%
    mutate(
      original_feature_id = feature_id,
      feature_id          = gene_symbol,
      id_space            = "symbol"
    ) %>%
    dplyr::select(-gene_symbol)
}


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

in_dir  <- file.path(per_study_dir, modality)
effect_files <- list.files(in_dir, pattern = "_effects\\.csv$", full.names = TRUE)
effect_files <- effect_files[!grepl("_normalized", effect_files)]

# Drop studies superseded by a finer-grained re-analysis of the same samples,
# which would otherwise be pooled twice. Each exclusion is logged.
superseded <- superseded_studies(
  file.path(REPO_ROOT, "conf", "analysis", "superseded_studies.csv"))
if (length(superseded) > 0) {
  accessions <- sub("_effects\\.csv$", "", basename(effect_files))
  drop <- accessions %in% superseded
  for (acc in accessions[drop]) {
    message(sprintf("Excluding %s: superseded (see conf/analysis/superseded_studies.csv)",
                    acc))
  }
  effect_files <- effect_files[!drop]
}

if (length(effect_files) == 0) {
  message("No effect files found in ", in_dir)
  quit(status = 0)
}

all_normalized <- list()

for (eff_path in effect_files) {
  accession <- sub("_effects\\.csv$", "", basename(eff_path))
  message(sprintf("\n--- %s ---", accession))

  df <- tryCatch(
    read.csv(eff_path, stringsAsFactors = FALSE),
    error = function(e) { message("  Read error: ", e$message); NULL }
  )
  if (is.null(df) || nrow(df) == 0) next
  # Some platforms use bare numeric probe ids (e.g. 6602010), which read.csv
  # types as double while every other study yields character. bind_rows then
  # refuses to combine the per-study frames, so pin the type on read.
  if ("feature_id" %in% colnames(df)) {
    df$feature_id <- as.character(df$feature_id)
    # Applied to the frame, not only to the lookup: the normalized table is
    # what 06 pools, so a study left in `gene-Dgkh` space would join nothing
    # even if the map behind it were right.
    stripped <- strip_gene_prefix(df$feature_id)
    if (!identical(stripped, df$feature_id)) {
      message("  Stripped the NCBI gene- prefix from ", length(stripped), " ids")
      df$feature_id <- stripped
    }
  }

  # First try GPL-based dispatch (Affymetrix / Illumina probe ID arrays)
  gpl <- detect_platform(accession, cache_dir)
  id_map <- NULL
  if (!is.na(gpl) && gpl %in% names(GPL_TO_ANNO_PKG)) {
    message("  GPL: ", gpl, " → ", GPL_TO_ANNO_PKG[[gpl]])
    id_map <- map_via_annotation_db(df$feature_id, GPL_TO_ANNO_PKG[[gpl]])
  }

  # Fall back to the GEO platform table, which is already cached for every
  # study. This recovers arrays whose feature ids are bare array-address
  # numbers (5617701, 10700001) rather than the string PROBEIDs a .db package
  # keys on, and arrays that have no Bioconductor package at all. Before this,
  # six studies and ~370,000 features sat in probe space in the published
  # results -- 56% of pooled features were probe ids, including 409 of the
  # 1,178 significant ones -- because the .db lookup was the only route.
  if (is.null(id_map)) {
    id_map <- map_via_gpl_table(df$feature_id, gpl, cache_dir)
  }

  # Fall back to ID-space heuristics (Ensembl transcripts, MTA-1.0)
  if (is.null(id_map)) {
    id_space <- detect_id_space(df$feature_id)
    message("  Detected ID space: ", id_space, " (GPL=", gpl, ")")
    id_map <- switch(id_space,
      rat_ensembl_transcript  = map_rat_ensembl_transcript(df$feature_id, interim_dir),
      mouse_mta10             = map_mouse_mta10(df$feature_id, interim_dir),
      mouse_ensembl_transcript = {
        map_ensembl_gene(df$feature_id, "mouse", interim_dir)
      },
      mouse_ensembl_gene      = map_ensembl_gene(df$feature_id, "mouse", interim_dir),
      rat_ensembl_gene        = map_ensembl_gene(df$feature_id, "rat",   interim_dir),
      human_ensembl_gene      = map_ensembl_gene(df$feature_id, "human", interim_dir),
      human_ensembl_transcript = map_ensembl_transcript(df$feature_id, "human",
                                                        interim_dir),
      entrez                  = map_entrez(df$feature_id, species_for(accession)),
      refseq                  = map_refseq(df$feature_id, species_for(accession)),
      {
        message("  No normalization needed or unsupported ID space — skipping")
        NULL
      }
    )
  }

  if (is.null(id_map)) {
    message("  Could not obtain ID map — keeping native IDs")
    df$original_feature_id <- df$feature_id
    df <- canonicalize_study_symbols(df, accession)$df
    all_normalized[[accession]] <- df
    # Write it anyway. 06 and 06g glob *_effects_normalized.csv and prefer it
    # over _effects.csv, so skipping the write does not make them fall back --
    # it leaves whatever normalized file the study had last time, which then
    # shadows the current effects. After the 2026-08-28 re-run that was 19
    # studies still serving May tables to the meta-analysis. Writing every run
    # also carries the symbol canonicalisation through to 06, which previously
    # reached only the combined file for these studies.
    out_path <- file.path(in_dir, paste0(accession, "_effects_normalized.csv"))
    write.csv(df, out_path, row.names = FALSE)
    message(sprintf("  %d rows saved (native ids) → %s", nrow(df), basename(out_path)))
    next
  }

  norm_df <- aggregate_to_gene(df, id_map)
  if (is.null(norm_df) || nrow(norm_df) == 0) {
    message("  No features mapped to gene symbols — skipping")
    next
  }
  norm_df <- canonicalize_study_symbols(norm_df, accession)$df

  out_path <- file.path(in_dir, paste0(accession, "_effects_normalized.csv"))
  write.csv(norm_df, out_path, row.names = FALSE)
  message(sprintf("  %d genes saved → %s", nrow(norm_df), basename(out_path)))
  all_normalized[[accession]] <- norm_df
}

# Combined output
if (length(all_normalized) > 0) {
  combined <- bind_rows(all_normalized)
  combined_path <- file.path(in_dir, "effects_normalized_all.csv")
  write.csv(combined, combined_path, row.names = FALSE)
  message(sprintf("\nCombined normalized file: %d rows from %d studies → %s",
                  nrow(combined), length(all_normalized), combined_path))
} else {
  message("\nNo studies successfully normalized.")
}

message("\nStep 05b complete.")
