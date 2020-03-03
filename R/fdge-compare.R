#' @rdname fdge
#' @export
#'
#' @section Comapring DGE Results:
#' We can compare two Ttest results.
#'
#' The filtering strategy in the interaction model dictates that the union
#' of all features found in `x` are `y` are used in the test.
#'
#' @param rerun When comparing two results, the features analyzed in each may
#'   differ, making comparisons between the two objects sparse, at times.
#'   When `rerun = TRUE` (default), the original linear models are rerun with
#'   their features set to `union(features(x), features(y))`.
#'
#' @examples
#' # Comparing two T-test results ----------------------------------------------
#' # Let's compare the tumor vs normal DGE results in CRC vs BLCA
#'
#' efds <- exampleFacileDataSet()
#' dge.crc <- filter_samples(efds, indication == "CRC") %>%
#'   flm_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.blca <- filter_samples(efds, indication == "BLCA") %>%
#'   flm_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.comp <- compare(dge.crc, dge.blca)
#' if (interactive()) {
#'   report(dge.comp)
#'   shine(dge.comp)
#' }
compare.FacileTtestAnalysisResult <- function(x, y, treat_lfc = NULL,
                                              rerun = TRUE, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  # TODO: assert_comparable(x, y, ...)
  # TODO: assert_comparable.FacileDgeAnalysisResult <- function(...)
  assert_class(x, "FacileTtestAnalysisResult")
  assert_class(y, "FacileTtestAnalysisResult")
  fds. <- assert_class(fds(x), "FacileDataStore")

  stopifnot(
    param(x, "assay_name") == param(y, "assay_name"),
    param(x, "method") == param(y, "method"))

  if (!is.null(treat_lfc)) {
    if (!test_number(treat_lfc, lower = 0)) {
      warnings <- c(
        warnings,
        "Illegal parameter passed to `treat_lfc`. It is being ignored")
      treat_lfc <- 0
    }
  } else {
    treat_lfc <- 0
  }

  # override downstream with "failed/incomplete/with_i_stats" version of class?
  clazz <- NULL
  classes <- c("FacileTtestComparisonAnalysisResult",
               "FacileTtestAnalysisResult",
               "FacileDgeAnalysisResult",
               "FacileComparisonAnalysis",
               "FacileAnalysisResult")

  out <- list(
    result = NULL,
    xystats = NULL,
    params = list(x = x, y = y, treat_lfc = treat_lfc),
    fds = fds.)

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, classes)
    return(out)
  })

  idge <- .interaction_fdge(x, y, treat_lfc = treat_lfc, rerun = rerun)
  if (is.null(idge)) {
    samples. <- bind_rows(samples(x), samples(y))
    samples. <- set_fds(samples., fds.)
    samples. <- distinct(samples., dataset, sample_id)
  } else {
    classes <- c("FacileTtestComparisonInteractionAnalysisResult", classes)
    samples. <- samples(idge[["result"]])
    if (idge[["rerun"]]) {
      # This is non-orthodox-facile:
      # 1. the x,y params may not be the precise versions that were sent in here
      out[["params"]][["x"]] <- x <- idge[["x"]]
      out[["params"]][["y"]] <- y <- idge[["y"]]
    }
  }

  xres <- tidy(x)
  yres <- tidy(y)

  jcols <- intersect(colnames(xres), colnames(yres))

  meta.cols <- c("feature_type", "feature_id", "symbol", "meta")
  drop.cols <- c("seqnames", "start", "end", "strand", "effective_length",
                 "source")
  stat.cols <- setdiff(colnames(xres),  c(meta.cols, drop.cols))

  meta.cols <- intersect(meta.cols, jcols)
  xystats <- full_join(
    select(xres, {{meta.cols}}, {{stat.cols}}),
    select(yres, feature_type, feature_id, {{stat.cols}}),
    by = c("feature_type", "feature_id"))

  if (!is.null(idge)) {
    ires <- tidy(idge[["result"]]) %>%
      select(feature_type, feature_id, {{stat.cols}})
    xystats <- left_join(xystats, ires, by = c("feature_type", "feature_id"))
    # put stats for interaction test up front, followed by *.x, *.y
    # xystats <- select(xystats, !!c(meta.cols, stat.cols), everything())
    xystats <- select(xystats, {{meta.cols}}, {{stat.cols}}, everything())
  }

  out[["result"]] <- idge[["result"]]
  out[["xystats"]] <- xystats
  out[["samples"]] <- samples.

  out
}

#' This `result()` of comparison will return the FacileTtestAnalysisResult
#' for the interaction, if it was run, otherwise NULL.
#'
#' This definition is redundant, since result.FacileAnalysisResult x returns
#' x[["result"]], but I just want to make this explicit for consumers of this
#' code (even myself).
#'
#' @noRd
#' @export
result.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  x[["result"]]
}

#' @noRd
#' @export
tidy.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  x[["xystats"]]
}

#' @noRd
#' @export
model.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  # model.NULL takes care of if/else testing NULL-ness of result()
  model(result(x))
}

#' @noRd
#' @export
samples.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  x[["samples"]]
}

#' @noRd
#' @export
viz.FacileTtestComparisonAnalysisResult <- function(x, max_padj = 0.1,
                                                    feature_id = NULL, ...) {
  hover <- c(
    # feature metadata
    "symbol", "feature_id", "meta",
    # padj from individual fdge tests
    "padj.x", "padj.y",
    # stats from interaction model, if it was run
    "logFC", "padj")
  d <- tidy(x)

  if (!is.null(feature_id)) {
    assert_multi_class(feature_id, c("character", "data.frame"))
    if (is.data.frame(feature_id)) {
      .feature_id <- feature_id[["feature_id"]]
    } else {
      .feature_id <- feature_id
    }
    d <- filter(d, feature_id %in% .feature_id)
  } else if (!is.null(x[["dge"]])) {
    d <- filter(d, padj.x <= max_padj | padj.y <= max_padj | padj <= max_padj)
  } else {
    d <- filter(d, padj.x <= max_padj | padj.y <= max_padj)
  }

  hover <- intersect(hover, colnames(d))
  fscatterplot(d, c("logFC.x", "logFC.y"), hover = hover, webgl = TRUE, ...)
}

#' Helper function to run an interaction model to generate statistics when
#' comparing two Ttest models.
#'
#' This function is only to be called from within compare.FacileTTestDGEResult.
#'
#' If we can't generate an interaction result, this will return NULL.
#' @noRd
.interaction_fdge <- function(x, y, treat_lfc = NULL, rerun = FALSE, ...) {
  # If these results aren't from the same FacileDataStore, get outta here
  # NOTE: This is not really a robust way to compare if two fds are the same
  if (name(fds(x)) != name(fds(y))) {
    return(NULL)
  }
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

  warning("Properly running the interaction model here is still ALPHA\n  ",
          "https://github.com/facileverse/FacileAnalysis/issues/19",
          immediate. = TRUE)

  xsamples <- samples(xmod)
  ysamples <- samples(ymod)

  ibatch <- unique(c(param(xmod, "batch"), param(ymod, "batch")))

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
    # bind_rows can not be overriden, so we have to ensure that fds is stored
    # back into the samples. facile_frame
    samples. <- set_fds(bind_rows(xsamples, ysamples), fds(x))
    samples. <- distinct(samples., dataset, sample_id, .keep_all = TRUE)
  }

  if (length(ibatch)) {
    # remove and add values to be thorough
    for (bcov in ibatch)  samples.[[bcov]] <- NULL
    samples. <- with_sample_covariates(samples., ibatch)
  }

  # The x and y samples that came in here should have no samples that have
  # NA values in the covariates used under their testing scenario, but it's
  # possible that now that they are rbind'd together, they do.
  #
  # I *believe* this should only happen he xysamples samples. data.frame can
  # have rows with NA values for covariates that come from the `batch`
  # covaraites used in either x or y. If this is the case, I guess we just have
  # to remove that covariate from the model(?)
  test.na <- unique(c(icovariate, ibatch))
  has.na <- sapply(samples.[, test.na], function(vals) any(is.na(vals)))
  has.na <- names(has.na)[has.na]
  if (length(has.na)) {
    msg <- paste("The following covariates have samples with NA values, and",
                 "therefore can't be used in the interaction model: ",
                 paste(has.na, collapse = ","))
    not.batch <- setdiff(has.na, ibatch)
    if (length(not.batch)) {
      msg <- glue(
        msg, "\n\n",
        "These covariates are not just the 'batch' covarites in the upstream ",
        "`fdge` results. Skipping the interaction model ...",
        paste(not.batch, collapse = ","))
      warning(msg)
      return(NULL)
    }
    msg <- glue(
      msg, "\n\n",
      "The covariates HAVE BEEN REMOVED in order to run the ",
      "interaction fdge")
    warning(msg)
    ibatch <- setdiff(ibatch, has.na)
  }

  xcontrast <- xmod[["contrast_string"]]
  ycontrast <- ymod[["contrast_string"]]
  if (!same.samples) {
    xcontrast <- .prefix_contrast(xcontrast, "xgrp.")
    ycontrast <- .prefix_contrast(ycontrast, "ygrp.")
  }
  contrast. <- glue("( {xcontrast} ) - ( {ycontrast} )")

  imodel <- flm_def(samples., icovariate, batch = ibatch,
                    contrast. = contrast.)
  genes. <- unique(c(xres[["feature_id"]], yres[["feature_id"]]))

  ires <- fdge(imodel, features = genes.,
               method = param(x, "method"),
               assay_name = param(x, "assay_name"),
               with_sample_weights = param(x, "with_sample_weights"),
               treat_lfc = treat_lfc)

  rerun <- rerun && !setequal(xres[["feature_id"]], genes.)
  if (rerun) {
    x <- fdge(xmod, features = genes., method = param(x, "method"),
              assay_name = param(x, "assay_name"),
              with_sample_weights = param(x, "with_sample_weights"))
    y <- fdge(ymod, features = genes., method = param(y, "method"),
              assay_name = param(y, "assay_name"),
              with_sample_weights = param(y, "with_sample_weights"))
  }

  list(result = ires, x = x, y = y, rerun = rerun)
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

#' @noRd
#' @export
print.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacileTtestComparisonAnalysisResult <- function(x, ...) {
  has.istats <- !is.null(result(x))
  status <- if (has.istats) "with interaction" else "no interaction"

  out <- paste(
    "===========================================================\n",
    sprintf("FacileTtestComparisonAnalysisResult (%s)\n", status),
    "-----------------------------------------------------------\n",
    sep = "")

  if (has.istats) {
    mdef <- model(x)
    des <- design(mdef)
    formula <- mdef[["design_formula"]]
    test <- mdef[["contrast_string"]]
    res <- tidy(x)
    ntested <- nrow(res)
    nsig <- sum(!is.na(res[["padj"]]) & res[["padj"]] < 0.10)
    out <- paste(
      out,
      glue("Significant Results (FDR < 0.1): ({nsig} / {ntested})"), "\n",
      "Formula: ", formula, "\n",
      "Tested: ", test, "\n",
      sep = "")
  } else {
    xform <- model(param(x, "x"))[["contrast_string"]]
    yform <- model(param(x, "y"))[["contrast_string"]]
    out <- paste(out, sprintf("(%s) - (%s)\n", xform, yform), sep = "")
  }

  paste(out,
        "===========================================================\n",
        sep = "")
}
