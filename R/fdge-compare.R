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
compare.FacileTtestDGEResult <- function(x, y, ...) {
  if (FALSE) {
    x <- dge.crc; y <- dge.blca
  }
  # TODO: assert_comparable(x, y, ...)
  # TODO: assert_comparable.FacileDGEResult <- function(...)
  assert_class(x, "FacileTtestDGEResult")
  assert_class(y, "FacileTtestDGEResult")
  stopifnot(
    param(x, "assay_name") == param(y, "assay_name"),
    param(x, "method") == param(y, "method"))

  fds. <- fds(x)
  xres <- result(x)
  yres <- result(y)

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
    ires <- select(result(idge), !!c(meta.cols, stat.cols))
    xystats <- left_join(xystats, ires, by = meta.cols)
  }

  out <- list(
    result = xystats,
    dge = idge)
  class(out) <- c("FacileTtestComparison",
                  "FacileAnalysisComparison",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
report.FacileTtestComparison <- function(x, ...) {
  if (!is.null(x[["dge"]])) {
    fscatterplot(result(x),
                 c("logFC.x", "logFC.y"),
                 hover = c("padj.x", "padj.y", "logFC", "padj"),
                 webgl = TRUE)
  } else {
    fscatterplot(result(x),
                 c("logFC.x", "logFC.y"),
                 hover = c("padj.x", "padj.y"),
                 webgl = TRUE)
  }
}

#' Helper function to run an interaction model to generate statistics when
#' comparing two Ttest models.
#'
#' This function is only to be called from within compare.FacileTTestDGEResult.
#'
#' If we can't generate an interaction result, this will return NULL.
.interaction_fdge <- function(x, y, ...) {
  xmod <- model(x)
  ymod <- model(y)
  xres <- result(x)
  yres <- result(y)

  covariate <- param(xmod, "covariate")
  if (covariate != param(ymod, "covariate")) {
    warning("Covariates used in test are not equal, no dge analysis performed",
            immediate. = TRUE)
    return(NULL)
  }

  icovariate <- ".grp."
  ifixed <- unique(c(param(xmod, "fixed"), param(ymod, "fixed")))
  xsamples <- samples(xmod)
  xadd <- setdiff(ifixed, colnames(xsamples))
  if (length(xadd)) {
    warning("Adding fixed covariates to x, these were not originally used: ",
            paste(xadd, collapse = ", "))
    xsamples <- with_sample_covariates(xsamples, .fds = fds(x))
  }
  xsamples[[icovariate]] <- paste0("xgrp.", xsamples[[covariate]])

  ysamples <- samples(ymod)
  yadd <- setdiff(ifixed, colnames(ysamples))
  if (length(yadd)) {
    warning("Adding fixed covariates to y, these were not originally used: ",
            paste(yadd, collapse = ", "))
    ysamples <- with_sample_covariates(ysamples, .fds = fds(y))
  }
  ysamples[[icovariate]] <- paste0("ygrp.", ysamples[[covariate]])

  samples. <- set_fds(bind_rows(xsamples, ysamples), fds(xmod))
  browser()
  contrast. <- glue(
    "( {xcon} ) - ( {ycon} )",
    xcon = .prefix_contrast(xmod[["contrast_string"]], "xgrp."),
    ycon = .prefix_contrast(ymod[["contrast_string"]], "ygrp."))

  imodel <- fdge_model_def(samples., icovariate, fixed = ifixed,
                           contrast. = contrast.)
  genes. <- unique(c(xres[["feature_id"]], yres[["feature_id"]]))

  fdge(imodel, filter = genes.,
       method = param(x, "method"),
       assay_name = param(x, "assay_name"))
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

