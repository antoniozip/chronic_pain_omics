#' Forest-plot drawing shared by step 06 and the manuscript figure.
#'
#' Split out of `pipeline/06_meta_analysis.R` because the two artefacts drawn
#' from these panels have opposite needs and were being served by one file:
#'
#'   - the manuscript figure shows a single feature and is placed in an 8.2 cm
#'     column, so its page must be the size of that one panel;
#'   - the supplementary set shows the top N features, whose panels differ in
#'     depth, so each of its pages is sized for its own panel.
#'
#' Serving both from one file gave the manuscript a page sized for the deepest
#' feature in the corpus: RCC2-AS1 has 5 contributing studies and arrived on a
#' 10x30 inch page, 85% of it blank. LaTeX then reported "Float too large for
#' page by 101 pt", deferred it, and a deferred float blocks every later float
#' of its class -- which is why every figure in the paper was landing on pages
#' 24-25 rather than beside the text that discusses it.
#'
#' No library dependencies beyond metafor, so this file can be sourced by a
#' regeneration script without pulling in the pipeline.

#' Page height in inches for a forest panel of `k` studies.
#'
#' The two summary rows (the header rule and the pooled diamond) and the axis
#' need room whatever `k` is, which is what the constant term covers; below
#' about four studies the panel is dominated by them, hence the floor.
forest_page_height <- function(k) {
  max(3.2, k * 0.35 + 2)
}


#' Draw one feature's random-effects forest panel.
#'
#' `feat_df` must already be filtered to usable standard errors: refitting the
#' rows `meta_one_feature()` used is what keeps the drawn diamond equal to the
#' pooled estimate reported for the feature. Returns TRUE if a panel was
#' drawn, FALSE if the feature had too few usable rows or the fit failed.
draw_forest_panel <- function(feat_df, feature_id, method = "REML",
                              slab_max_chars = 60) {
  if (nrow(feat_df) < 2) return(FALSE)
  ok <- TRUE
  tryCatch({
    res <- metafor::rma(yi = effect_size, sei = se, data = feat_df,
                        method = method, slab = feat_df$study_id)
    metafor::forest(res,
                    main = substr(feature_id, 1, slab_max_chars),
                    xlab = unique(feat_df$effect_unit)[1])
  }, error = function(e) ok <<- FALSE)
  ok
}


#' Combine single-page PDFs into one file, each page keeping its own size.
#'
#' A `pdf()` device fixes one page size for everything it writes, so a set
#' whose panels differ in depth has to be drawn a page at a time and joined
#' afterwards. Returns TRUE if the parts were joined, FALSE if no tool is
#' available or the join failed, which lets the caller fall back rather than
#' leave no figure at all.
merge_pdf_pages <- function(parts, out_path) {
  pdfunite <- Sys.which("pdfunite")
  if (nzchar(pdfunite)) {
    ok <- system2(pdfunite, shQuote(c(parts, out_path)),
                  stdout = FALSE, stderr = FALSE)
    if (identical(ok, 0L)) return(TRUE)
  }
  ghostscript <- Sys.which("gs")
  if (nzchar(ghostscript)) {
    # -dAutoRotatePages=/None, or a landscape panel is turned upright and its
    # page size silently swapped.
    ok <- system2(ghostscript,
                  c("-dNOPAUSE", "-dBATCH", "-dQUIET", "-sDEVICE=pdfwrite",
                    "-dAutoRotatePages=/None",
                    paste0("-sOutputFile=", shQuote(out_path)), shQuote(parts)),
                  stdout = FALSE, stderr = FALSE)
    if (identical(ok, 0L)) return(TRUE)
  }
  FALSE
}


#' Write forest panels to a PDF, one page per feature.
#'
#' `features` is drawn in the order given, each on a page sized to its own
#' panel: the pages of a set differ in height, because the features differ in
#' how many studies measured them. That needs one device per page and a join,
#' since a single `pdf()` device would impose the deepest panel's page on all
#' of them -- which is how the 50-page transcriptomic set came to be 50 pages
#' of 34 inches for a median panel of 4, `forest()` stretching five studies
#' down a page sized for ninety.
#'
#' Without a tool to join them the whole set goes on one page size, as before,
#' with a warning: a stretched figure beats no figure.
write_forest_pdf <- function(out_path, effects_df, features, method = "REML",
                             slab_max_chars = 60, width = 10) {
  subset_df <- effects_df[effects_df$feature_id %in% features, , drop = FALSE]
  if (nrow(subset_df) == 0) return(invisible(NULL))
  features <- features[features %in% subset_df$feature_id]

  panel_df <- function(feat) {
    usable_se(subset_df[subset_df$feature_id == feat, , drop = FALSE])
  }
  draw_page <- function(feat, path) {
    feat_df <- panel_df(feat)
    grDevices::pdf(path, width = width,
                   height = forest_page_height(nrow(feat_df)))
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_forest_panel(feat_df, feat, method = method,
                      slab_max_chars = slab_max_chars)
  }

  if (length(features) == 1) {
    draw_page(features, out_path)
    return(invisible(out_path))
  }

  parts_dir <- tempfile("forest_pages")
  dir.create(parts_dir)
  on.exit(unlink(parts_dir, recursive = TRUE), add = TRUE)
  parts <- character(0)
  for (i in seq_along(features)) {
    part <- file.path(parts_dir, sprintf("page_%03d.pdf", i))
    # A feature whose panel could not be drawn contributes no page, rather
    # than a pageless part the join would refuse.
    if (isTRUE(draw_page(features[i], part))) parts <- c(parts, part)
  }
  if (length(parts) > 0 && merge_pdf_pages(parts, out_path)) {
    return(invisible(out_path))
  }

  warning("no pdfunite or gs to join the forest pages: writing ", out_path,
          " on one page size, sized for its deepest panel", call. = FALSE)
  max_studies <- max(vapply(features, function(f) nrow(panel_df(f)), integer(1)))
  grDevices::pdf(out_path, width = width,
                 height = forest_page_height(max_studies))
  on.exit(grDevices::dev.off(), add = TRUE)
  for (feat in features) {
    draw_forest_panel(panel_df(feat), feat, method = method,
                      slab_max_chars = slab_max_chars)
  }
  invisible(out_path)
}
