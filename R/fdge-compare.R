#' @section Comapring DGE Results:
#' We can compare two Ttest results
#'
#' @rdname fdge
#' @export
#' @examples
#' # Comparing two T-test results ----------------------------------------------
#' # Let's compare the tumor vs normal DGE results in CRC vs BLCA
#'
#' efds <- exampleFacileDataSet()
#' dge.crc <- filter_samples(efds, indication == "CRC") %>%
#'   fdge_model_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.blca <- filter_samples(efds, indication == "BLCA") %>%
#'   fdge_model_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.comp <- compare(dge.crc, dge.blca)
#' if (interactive()) {
#'   report(dge.comp)
#'   shine(dge.comp)
#' }
compare.FacileTtestAnalysisResult <- function(x, y, ...) {
  # TODO: assert_comparable(x, y, ...)
  # TODO: assert_comparable.FacileDgeAnalysisResult <- function(...)
  assert_class(x, "FacileTtestAnalysisResult")
  assert_class(y, "FacileTtestAnalysisResult")
  stopifnot(
    param(x, "assay_name") == param(y, "assay_name"),
    param(x, "method") == param(y, "method"))

  fds. <- fds(x)
  xres <- tidy(x)
  yres <- tidy(y)

  idge <- .interaction_fdge(x, y)

  meta.cols <- c("feature_type", "feature_id", "symbol", "meta")
  drop.cols <- c("seqnames", "start", "end", "strand", "effective_length",
                 "source")
  stat.cols <- setdiff(colnames(xres),  c(meta.cols, drop.cols))
  xystats <- full_join(
    select(xres, !!c(meta.cols, stat.cols)),
    select(yres, !!c(meta.cols, stat.cols)),
    by = meta.cols)

  if (!is.null(idge)) {
    ires <- select(tidy(idge), !!c(meta.cols, stat.cols))
    xystats <- left_join(xystats, ires, by = meta.cols)
    # put stats for interaction test up front, followed by *.x, *.y
    xystats <- select(xystats, !!c(meta.cols, stat.cols), everything())
  }

  out <- list(
    result = xystats,
    dge = idge)
  class(out) <- c("FacileTtestComparisonAnalysisResult",
                  "FacileTtestAnalysisResult",
                  "FacileAnalysisComparisonComparison",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
report.FacileTtestComparisonAnalysisResult <- function(x, max_padj = 0.1, ...) {
  hover <- c(
    # feature metadata
    "symbol", "feature_id", "meta",
    # padj from individual fdge tests
    "padj.x", "padj.y",
    # stats from interaction model, if it was run
    "logFC", "padj")
  d <- tidy(x)
  if (!is.null(x[["dge"]])) {
    d <- filter(d, padj.x <= max_padj | padj.y <= max_padj | padj <= max_padj)
  } else {
    d <- filter(d, padj.x <= max_padj | padj.y <= max_padj)
  }
  hover <- intersect(hover, colnames(d))
  fscatterplot(d, c("logFC.x", "logFC.y"), hover = hover, webgl = TRUE)
}

#' Helper function to run an interaction model to generate statistics when
#' comparing two Ttest models.
#'
#' This function is only to be called from within compare.FacileTTestDGEResult.
#'
#' If we can't generate an interaction result, this will return NULL.
#' @noRd
.interaction_fdge <- function(x, y, ...) {
  xmod <- model(x)
  ymod <- model(y)
  xres <- tidy(x)
  yres <- tidy(y)

  covariate <- param(xmod, "covariate")
  if (covariate != param(ymod, "covariate")) {
    warning("Covariates used in test are not equal, no dge analysis performed",
            immediate. = TRUE)
    return(NULL)
  }

  warning("Properly running the interaction model here is still ALPHA",
          immediate. = TRUE)

  xsamples <- samples(xmod)
  ysamples <- samples(ymod)

  ifixed <- unique(c(param(xmod, "fixed"), param(ymod, "fixed")))

  # If the models were defined on separate samples, then we have to take
  # extra care, since if the samples were split differently, then the same
  # covariate values can identify the same samples across models.
  # For example, if we did a tumor vs normal comparison in x of CRC and the
  # same comparison in y from BLCA, then we need to be more careful with
  # contrast for (crc tumor / crc normal) / (blca tumor / blca normal)
  same.samples <- setequal(
    paste(xsamples$dataset, xsamples$sample_id),
    paste(ysamples$dataset, ysamples$sample_id))

  if (same.samples) {
    samples. <- xsamples
    icovariate <- covariate
  } else {
    icovariate <- ".grp."
    xsamples[[icovariate]] <- paste0("xgrp.", xsamples[[covariate]])
    ysamples[[icovariate]] <- paste0("ygrp.", ysamples[[covariate]])
    samples. <- bind_rows(xsamples, ysamples)
    samples. <- distinct(samples., dataset, sample_id, .keep_all = TRUE)
  }

  samples. <- set_fds(samples., fds(x))

  if (length(ifixed)) {
    # remove and add values to be thorough
    for (fixcov in ifixed)  samples.[[fixcov]] <- NULL
    samples. <- with_sample_covariates(samples., ifixed)
  }

  xcontrast <- xmod[["contrast_string"]]
  ycontrast <- ymod[["contrast_string"]]
  if (!same.samples) {
    xcontrast <- .prefix_contrast(xcontrast, "xgrp.")
    ycontrast <- .prefix_contrast(ycontrast, "ygrp.")
  }
  contrast. <- glue("( {xcontrast} ) - ( {ycontrast} )")

  imodel <- fdge_model_def(samples., icovariate, fixed = ifixed,
                           contrast. = contrast.)
  genes. <- unique(c(xres[["feature_id"]], yres[["feature_id"]]))

  fdge(imodel, filter = genes.,
       method = param(x, "method"),
       assay_name = param(x, "assay_name"),
       with_sample_weights = param(x, "with_sample_weights"))
}

#' @noRd
#' @importFrom stringr str_split
.prefix_contrast <- function(eqn, prefix = "x.") {
  assert_string(eqn)
  assert_string(prefix)
  splitted <- str_split(eqn, " +")[[1L]]
  out <- lapply(splitted, function(x) {
    if (x == make.names(x)) paste0(prefix, x) else x
  })
  paste(out, collapse = " ")
}

