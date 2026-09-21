#' Forest-plot drawing shared by step 06 and the manuscript figure.
#'
#' Split out of `pipeline/06_meta_analysis.R` because the two artefacts drawn
#' from these panels have opposite needs and were being served by one file:
#'
#'   - the manuscript figure shows a single feature and is placed in an 8.2 cm
#'     column, so its page must be the size of that one panel;
#'   - the supplementary set shows the top N features, and a `pdf()` device
#'     fixes one page size for the whole file, so its pages are sized for the
#'     deepest panel in the set.
#'
#' Serving both from the multi-page file gave the manuscript a page sized for
#' the deepest feature in the corpus: RCC2-AS1 has 3 contributing studies and
#' arrived on a 10x30 inch page, 85% of it blank. LaTeX then reported "Float
#' too large for page by 101 pt", deferred it, and a deferred float blocks
#' every later float of its class -- which is why every figure in the paper
#' was landing on pages 24-25 rather than beside the text that discusses it.
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


#' Write forest panels to a PDF, one page per feature.
#'
#' `features` is drawn in the order given. The page size is fixed for the file,
#' so it is taken from the deepest panel in the set: pass one feature to get a
#' page that fits it exactly, which is what the manuscript figure needs.
write_forest_pdf <- function(out_path, effects_df, features, method = "REML",
                             slab_max_chars = 60, width = 10) {
  subset_df <- effects_df[effects_df$feature_id %in% features, , drop = FALSE]
  if (nrow(subset_df) == 0) return(invisible(NULL))

  max_studies <- max(table(subset_df$feature_id))
  grDevices::pdf(out_path, width = width,
                 height = forest_page_height(max_studies))
  on.exit(grDevices::dev.off(), add = TRUE)

  for (feat in features) {
    feat_df <- usable_se(subset_df[subset_df$feature_id == feat, , drop = FALSE])
    draw_forest_panel(feat_df, feat, method = method,
                      slab_max_chars = slab_max_chars)
  }
  invisible(out_path)
}
