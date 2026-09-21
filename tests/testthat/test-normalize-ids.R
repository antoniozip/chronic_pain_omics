library(testthat)

# 05b_normalize_ids.R cannot be source()d here: its top-level block loads
# GEOquery and the Bioconductor annotation packages. As in
# test-transcriptomics-da.R, lift individual functions out of the source.

REPO_ROOT <- normalizePath(file.path(test_path(), "..", ".."))
SRC <- readLines(file.path(REPO_ROOT, "pipeline", "05b_normalize_ids.R"), warn = FALSE)
CODE <- SRC[!grepl("^\\s*#", SRC)]

lift <- function(names) {
  ends <- grep("^\\}|^\\)", SRC)
  out <- character()
  for (n in names) {
    at <- grep(paste0("^", n, " *<-"), SRC)
    stopifnot(length(at) == 1)
    # A one-line definition ends on its own line; a multi-line one ends at the
    # next line starting with a closing brace or paren.
    if (grepl("^\\s*$|[\\)\\}\"]\\s*$", SRC[at]) && !grepl("\\($", trimws(SRC[at]))) {
      out <- c(out, SRC[at])
    } else {
      close_at <- ends[ends > at]
      out <- c(out, SRC[at:(if (length(close_at)) min(close_at) else at)])
    }
  }
  eval(parse(text = paste(out, collapse = "\n")), envir = parent.frame())
}


# ---------------------------------------------------------------------------
# detect_id_space() reads a threshold pair, a share helper and a pattern table,
# so lifting the function alone leaves those undefined.
ID_SPACE_FUNCTIONS <- c("ID_SPACE_MIN_SHARE", "ID_SPACE_UNANIMOUS",
                        "id_share", "ID_SPACE_PATTERNS", "detect_id_space")

# detect_id_space
#
# Five studies (137,154 rows) reached the meta-analysis keyed by something
# other than a gene symbol, because this function had branches only for mouse
# and rat Ensembl plus MTA-1.0. Everything else returned "unknown", which
# routes to "keeping native IDs" -- indistinguishable in the output from a
# study whose ids were already symbols.
# ---------------------------------------------------------------------------

test_that("detect_id_space recognises human Ensembl genes", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("ENSG00000000003", "ENSG00000000005")),
               "human_ensembl_gene")
})

test_that("detect_id_space recognises versioned human Ensembl genes", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("ENSG00000012779.11", "ENSG00000138600.10")),
               "human_ensembl_gene")
})

# GSE153739 and GSE153740 are 100% human Ensembl *transcript* accessions.
# detect_id_space() had ENSMUST and ENSRNOT branches but no ENST one, so both
# fell through to "unknown", kept their native ids, and entered the pool as
# 215,320 features that no other study can ever join -- 58% of the whole-cohort
# feature count, contributing nothing to the BH family.
test_that("detect_id_space recognises human Ensembl transcripts", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("ENST00000369152", "ENST00000551855")),
               "human_ensembl_transcript")
  expect_equal(detect_id_space("ENST00000369152.7"), "human_ensembl_transcript")
})

test_that("detect_id_space does not read rodent transcripts as human ones", {
  lift(ID_SPACE_FUNCTIONS)
  # ENSMUST/ENSRNOT begin ENSM/ENSR, so the ENST branch cannot capture them --
  # but the branch order is what guarantees it, and order is easy to lose.
  expect_equal(detect_id_space("ENSMUST00000000001"), "mouse_ensembl_transcript")
  expect_equal(detect_id_space("ENSRNOT00000000001"), "rat_ensembl_transcript")
})

test_that("map_ensembl_transcript reads the interim cache", {
  lift("map_ensembl_transcript")
  dir <- withr::local_tempdir()
  write.csv(
    data.frame(ensembl_transcript_id = c("ENST00000369152", "ENST00000551855"),
               external_gene_name    = c("BOLA1", "CPNE8")),
    file.path(dir, "human_enst2symbol.csv"), row.names = FALSE)

  out <- map_ensembl_transcript(c("ENST00000369152", "ENST00000551855"),
                                "human", dir)
  expect_equal(unname(out[["ENST00000369152"]]), "BOLA1")
  expect_equal(unname(out[["ENST00000551855"]]), "CPNE8")
})

test_that("map_ensembl_transcript strips the version suffix on lookup", {
  lift("map_ensembl_transcript")
  dir <- withr::local_tempdir()
  write.csv(
    data.frame(ensembl_transcript_id = "ENST00000369152",
               external_gene_name    = "BOLA1"),
    file.path(dir, "human_enst2symbol.csv"), row.names = FALSE)

  out <- map_ensembl_transcript("ENST00000369152.7", "human", dir)
  expect_equal(unname(out[["ENST00000369152.7"]]), "BOLA1")
})

test_that("map_ensembl_transcript declines a cache that covers almost nothing", {
  lift("map_ensembl_transcript")
  dir <- withr::local_tempdir()
  write.csv(
    data.frame(ensembl_transcript_id = "ENST00000000001",
               external_gene_name    = "SOMEGENE"),
    file.path(dir, "human_enst2symbol.csv"), row.names = FALSE)

  # Nine of ten ids absent from the cache: keeping the accessions is right,
  # writing them out as if they were symbols is the defect this guards.
  ids <- sprintf("ENST0000000%04d", 1:10)
  out <- map_ensembl_transcript(ids, "human", dir, min_match = 0.3)
  # Falls through to the org package, which does not know these ids either.
  expect_true(is.null(out) || sum(!is.na(out)) <= 1)
})

test_that("detect_id_space recognises bare Entrez ids", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("100034253", "100036765", "100049583")), "entrez")
})

test_that("detect_id_space keeps its existing branches", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space("ENSMUSG00000000001"), "mouse_ensembl_gene")
  expect_equal(detect_id_space("ENSRNOG00000000001"), "rat_ensembl_gene")
  expect_equal(detect_id_space("ENSRNOT00000000001"), "rat_ensembl_transcript")
  expect_equal(detect_id_space("ENSMUST00000000001"), "mouse_ensembl_transcript")
  expect_equal(detect_id_space(c("Actb", "Gapdh")), "unknown")
})

test_that("detect_id_space does not call mixed symbols Entrez", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("100034253", "Actb")), "unknown")
})

test_that("a handful of ENSG contaminants does not make a study Ensembl-keyed", {
  # GSE221921 is 99.8% HGNC symbols (TSPAN6, DPM1) with 42 ENSG-prefixed rows
  # out of 21,914 -- and those are not even valid human Ensembl ids
  # (ENSG10010139717.1). `any()` flipped the whole study into the Ensembl
  # branch, where aggregate_to_gene() would drop the 99.8% that cannot map.
  lift(ID_SPACE_FUNCTIONS)
  ids <- c(rep("TSPAN6", 500), "ENSG10010139717.1")
  expect_equal(detect_id_space(ids), "unknown")
})

test_that("a genuinely Ensembl-keyed human study is still detected", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(paste0("ENSG", sprintf("%011d", 1:200))),
               "human_ensembl_gene")
})

test_that("detect_id_space tolerates a few non-Entrez stragglers", {
  # GSE117526 carries 6 RefSeq accessions among 14,346 Entrez ids. Requiring
  # every id to be numeric let 0.04% of the study veto the mapping for the
  # other 99.96%, and it stayed in Entrez space.
  lift(ID_SPACE_FUNCTIONS)
  ids <- c(as.character(100034253:100034352), "NM_001277450")
  expect_equal(detect_id_space(ids), "entrez")
})


# ---------------------------------------------------------------------------
# ensembl_dataset_for
#
# map_ensembl_gene() resolved species as `if mouse then mouse else rat`, so a
# human study was silently queried against the rat mart. That returns almost
# nothing and is indistinguishable from a poorly annotated study.
# ---------------------------------------------------------------------------

test_that("ensembl_dataset_for resolves the three species", {
  lift("ensembl_dataset_for")
  expect_equal(ensembl_dataset_for("Mus musculus")[["dataset"]], "mmusculus_gene_ensembl")
  expect_equal(ensembl_dataset_for("Rattus norvegicus")[["dataset"]], "rnorvegicus_gene_ensembl")
  expect_equal(ensembl_dataset_for("Homo sapiens")[["dataset"]], "hsapiens_gene_ensembl")
})

test_that("ensembl_dataset_for declines an unknown species instead of defaulting", {
  lift("ensembl_dataset_for")
  expect_null(ensembl_dataset_for("Danio rerio"))
  expect_null(ensembl_dataset_for(NA))
  expect_null(ensembl_dataset_for(""))
})

test_that("ensembl_dataset_for declines a multi-species study", {
  # GSE236754 is recorded as Mus musculus + Rattus norvegicus. Picking either
  # one silently would map half the features to the wrong organism.
  lift("ensembl_dataset_for")
  expect_null(ensembl_dataset_for(c("Mus musculus", "Rattus norvegicus")))
})


# ---------------------------------------------------------------------------
# map_via_gpl_table: ORF carries the symbol on Arraystar platforms
# ---------------------------------------------------------------------------

test_that("map_via_gpl_table reads the ORF column as a symbol", {
  lift(c("ensure_gpl_soft", "map_via_gpl_table"))
  dir <- withr::local_tempdir()
  writeLines(c("^PLATFORM = GPLTEST", "!platform_table_begin",
               "ID\tTRANSCRIPT_TYPE\tACC\tORF",
               "PROBE1\tmRNA\tENSMUST01\tRarg",
               "PROBE2\tmRNA\tENSMUST02\tActb",
               "!platform_table_end"),
             file.path(dir, "GPLTEST.soft.gz"))
  m <- map_via_gpl_table(c("PROBE1", "PROBE2"), "GPLTEST", dir)
  expect_equal(unname(m[c("PROBE1", "PROBE2")]), c("Rarg", "Actb"))
})

test_that("map_via_gpl_table still prefers an explicit symbol column over ORF", {
  lift(c("ensure_gpl_soft", "map_via_gpl_table"))
  dir <- withr::local_tempdir()
  writeLines(c("!platform_table_begin",
               "ID\tORF\tGeneSymbol",
               "P1\twrong\tRight1",
               "!platform_table_end"),
             file.path(dir, "GPLTEST.soft.gz"))
  expect_equal(unname(map_via_gpl_table("P1", "GPLTEST", dir)["P1"]), "Right1")
})


# ---------------------------------------------------------------------------
# map_ensembl_via_orgdb
#
# biomaRt::useMart() returns HTTP 404 on this workstation, so the network being
# reachable does not make the Ensembl route work. The org.*.eg.db packages
# carry the same ENSEMBL -> SYMBOL mapping offline; this is tried first, and
# biomaRt remains the fallback for anything they miss.
# ---------------------------------------------------------------------------

test_that("map_ensembl_via_orgdb maps human Ensembl genes offline", {
  skip_if_not_installed("org.Hs.eg.db")
  lift("map_ensembl_via_orgdb")
  m <- map_ensembl_via_orgdb(c("ENSG00000141510", "ENSG00000012048"), "human")
  expect_equal(unname(m[["ENSG00000141510"]]), "TP53")
  expect_equal(unname(m[["ENSG00000012048"]]), "BRCA1")
})

test_that("map_ensembl_via_orgdb strips the version suffix", {
  skip_if_not_installed("org.Hs.eg.db")
  lift("map_ensembl_via_orgdb")
  m <- map_ensembl_via_orgdb("ENSG00000141510.17", "human")
  expect_equal(unname(m[["ENSG00000141510.17"]]), "TP53")
})

test_that("map_ensembl_via_orgdb declines when nothing resolves", {
  skip_if_not_installed("org.Hs.eg.db")
  lift("map_ensembl_via_orgdb")
  expect_null(map_ensembl_via_orgdb(paste0("ENSG9999999", 1:20), "human"))
})


# ---------------------------------------------------------------------------
# Every processed study must write its normalized table.
#
# 05b wrote _effects_normalized.csv only for studies it could map, and used
# `next` for the rest. 06 and 06g glob *_effects_normalized.csv and prefer it
# over _effects.csv, so a study that stops mapping keeps whatever normalized
# file it had -- 19 studies were shadowed by May/August-12 files after today's
# 05 re-run, and the fresh effects never reached the meta-analysis. It also
# meant canonicalisation (March7 -> Marchf7) reached only the combined file.
# ---------------------------------------------------------------------------

test_that("the keep-native-ids branch writes a normalized table", {
  at <- grep("Could not obtain ID map", CODE)
  expect_length(at, 1)
  nxt <- grep("^\\s*next\\s*$", CODE)
  branch_end <- min(nxt[nxt > at])
  writes <- grep("write.csv", CODE)
  expect_true(any(writes > at & writes < branch_end),
              info = "keep-native-ids must write before `next`, or 06 reads a stale file")
})


# ---------------------------------------------------------------------------
# Two id spaces the 2026-08-29 amendment run met and 05b could not resolve.
# Both would have survived as native ids -- indistinguishable in the output
# from a study whose ids were already symbols -- and pooled as single-study
# features that can never join another study. See
# plan/2026-08-29-amendment-da-run.md.
# ---------------------------------------------------------------------------

test_that("strip_gene_prefix removes the NCBI GFF gene- prefix", {
  lift("strip_gene_prefix")
  expect_equal(strip_gene_prefix(c("gene-Dgkh", "gene-Trpv1")), c("Dgkh", "Trpv1"))
})

test_that("strip_gene_prefix leaves a table that does not carry the prefix", {
  lift("strip_gene_prefix")
  ids <- c("Dgkh", "Trpv1", "gene-Scn9a")
  expect_equal(strip_gene_prefix(ids), ids)
})

test_that("strip_gene_prefix leaves rna- and cds- features alone", {
  # Stripping those yields a transcript or CDS accession, not a gene symbol,
  # so they must fall through to the normal detection rather than be renamed.
  lift("strip_gene_prefix")
  ids <- rep(c("rna-XM_006238698.4", "cds-XP_006238760.1"), 5)
  expect_equal(strip_gene_prefix(ids), ids)
})

test_that("detect_id_space recognises versioned RefSeq transcripts", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("NM_000015.2", "NM_001042.3", "XM_005248378.1")),
               "refseq")
})

test_that("detect_id_space recognises unversioned RefSeq transcripts", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("NM_000015", "NR_002196", "NM_001042")), "refseq")
})

test_that("detect_id_space accepts an underscore version separator", {
  # GSE193928 deposits NM_022892_2 alongside NM_004409; requiring a dot put it
  # at 88% against the 0.95 threshold, so a wholly RefSeq study read as unknown.
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("NM_022892_2", "NM_004409", "NR_046768")), "refseq")
})

test_that("detect_id_space does not call a symbol table RefSeq", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c("TSPAN6", "DPM1", "SCYL3", "NM_000015.2")), "unknown")
})

test_that("the RefSeq route is dispatched from the id-space switch", {
  at <- grep("refseq +=", CODE)
  expect_gte(length(at), 1)
  expect_true(any(grepl("map_refseq\\(", CODE)))
})

test_that("map_refseq declines when the species cannot be resolved", {
  lift(c("ensembl_dataset_for", "map_refseq"))
  expect_null(map_refseq(c("NM_000015.2"), c("Mus musculus", "Rattus norvegicus")))
})


# ---------------------------------------------------------------------------
# The platform table has to be fetched, not assumed present.
#
# map_via_gpl_table() bailed out when {gpl}.soft.gz was absent, and nothing in
# 05b ever downloaded one: detect_platform() only reads the GPL accession out
# of the series matrix. Platforms an earlier run happened to fetch worked;
# the thirteen the 2026-08-28 amendment introduced did not, and fifteen studies
# -- 900k rows on Affymetrix, Illumina, Agilent and Arraystar arrays -- were
# written out in probe space, which is indistinguishable from a study whose ids
# were already symbols.
# ---------------------------------------------------------------------------

test_that("ensure_gpl_soft returns a cached table without fetching", {
  lift("ensure_gpl_soft")
  dir <- withr::local_tempdir()
  path <- file.path(dir, "GPL570.soft.gz")
  writeLines("x", path)
  expect_equal(ensure_gpl_soft("GPL570", dir), path)
})

test_that("ensure_gpl_soft accepts the uncompressed and .annot forms", {
  lift("ensure_gpl_soft")
  dir <- withr::local_tempdir()
  writeLines("x", file.path(dir, "GPL96.soft"))
  expect_equal(ensure_gpl_soft("GPL96", dir), file.path(dir, "GPL96.soft"))
  dir2 <- withr::local_tempdir()
  writeLines("x", file.path(dir2, "GPL97.annot.gz"))
  expect_equal(ensure_gpl_soft("GPL97", dir2), file.path(dir2, "GPL97.annot.gz"))
})

test_that("ensure_gpl_soft declines a missing platform accession", {
  lift("ensure_gpl_soft")
  dir <- withr::local_tempdir()
  expect_true(is.na(ensure_gpl_soft(NA_character_, dir)))
  expect_true(is.na(ensure_gpl_soft(NULL, dir)))
})

test_that("map_via_gpl_table asks for the table before giving up on it", {
  at_fetch <- grep("path <- ensure_gpl_soft(gpl, cache_dir)", CODE, fixed = TRUE)
  expect_length(at_fetch, 1)
  # and no bare file.exists bail-out survives ahead of it
  at_fn <- grep("^map_via_gpl_table <- function", CODE)
  expect_gt(at_fetch, at_fn)
  expect_lt(at_fetch, at_fn + 5)
})


# A stray minority routing the whole study is the mistake `any()` makes, and it
# is silent: the study emerges with the right columns and almost no rows. The
# single-cell arm's two rat studies hit it -- GSE198608 is 0.72% ENSRNOG and
# went from 21,170 features to 2; GSE289659 is 13.7% and went from 18,142 to
# 378. Both are otherwise rat gene symbols.
test_that("a small ENSRNOG minority does not make a symbol study rat-Ensembl", {
  lift(ID_SPACE_FUNCTIONS)
  ids <- c(rep("Atf3", 993), rep("ENSRNOG00000000891", 7))
  expect_equal(detect_id_space(ids), "unknown")
  ids <- c(rep("Gap43", 863), rep("ENSRNOG00000000053", 137))
  expect_equal(detect_id_space(ids), "unknown")
})

test_that("a small ENSMUSG minority does not make a symbol study mouse-Ensembl", {
  lift(ID_SPACE_FUNCTIONS)
  expect_equal(detect_id_space(c(rep("Actb", 95), rep("ENSMUSG00000000001", 5))),
               "unknown")
})

# The other direction: two studies in the corpus are legitimately just over 90%
# one space, and a 0.95 threshold would send both to `unknown` and leave them
# keyed on probesets and Ensembl ids. The margin is wide -- no study anywhere in
# the corpus sits between 20% and 90% of any space.
test_that("a study that is 90% MTA-1.0 probesets is still routed there", {
  lift(ID_SPACE_FUNCTIONS)
  ids <- c(rep("TC0100006437.mm.1", 905), rep("Actb", 95))   # GSE147216
  expect_equal(detect_id_space(ids), "mouse_mta10")
})

test_that("a study that is 94% rat Ensembl is still routed there", {
  lift(ID_SPACE_FUNCTIONS)
  ids <- c(rep("ENSRNOG00000000053", 943), rep("Actb", 57))  # GSE245768
  expect_equal(detect_id_space(ids), "rat_ensembl_gene")
})

test_that("a unanimous space is routed without the mixed-table warning", {
  lift(ID_SPACE_FUNCTIONS)
  expect_silent(space <- detect_id_space(rep("ENSMUSG00000000001", 100)))
  expect_equal(space, "mouse_ensembl_gene")
})

test_that("a non-unanimous majority says how much it will fail to map", {
  lift(ID_SPACE_FUNCTIONS)
  expect_message(detect_id_space(c(rep("ENSRNOG00000000053", 90),
                                   rep("Actb", 10))),
                 "mixed table")
})

test_that("every id space is tested by a share, never by any() or all()", {
  # The regression this guards: both `any()` and `all()` have shipped here, in
  # opposite directions, and neither raises.
  at <- grep("^detect_id_space <- function", SRC)
  close_at <- grep("^\\}", SRC); close_at <- close_at[close_at > at][1]
  body <- SRC[at:close_at]
  expect_false(any(grepl("\\bany\\(grepl", body)))
  expect_false(any(grepl("\\ball\\(grepl", body)))
  expect_true(any(grepl("ID_SPACE_MIN_SHARE", body)))
})
