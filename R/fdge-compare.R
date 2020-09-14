#' @include fdge.R
NULL

#' @rdname fdge
#' @export
#'
#' @section Comparing DGE Results:
#' It is often useful to compare the results of two t-tests, and for many
#' experimental designs, this can be a an intuitive way to perform test an
#' interaction effect.
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
#' efds <- FacileData::exampleFacileDataSet()
#' dge.crc <- FacileData::filter_samples(efds, indication == "CRC") %>%
#'   flm_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.blca <- FacileData::filter_samples(efds, indication == "BLCA") %>%
#'   flm_def("sample_type", "tumor", "normal", "sex") %>%
#'   fdge()
#' dge.comp <- compare(dge.crc, dge.blca)
#'
#' if (interactive()) {
#'   viz(dge.comp, xlabel = "logFC(CRC)", ylabel = "logFC(BLCA)",
#'       highlight = head(tidy(signature(dge.comp)), 5))
#'   report(dge.comp)
#'   shine(dge.comp)
#' }
#'
#' # Static visualization generates the main "4-way" plot, as well as the
#' # facets for each category.
#' sviz <- viz(dge.comp, static = TRUE, labels = c(x = "CRC", y = "BLCA"),
#'             subtitle = "Tumor vs normal comparisons across indications")
#' if (requireNamespace("patchwork")) {
#'   patchwork::wrap_plots(
#'     sviz$plot + ggplot2::theme(legend.position = "bottom"),
#'     sviz$plot_facets + ggplot2::theme(legend.position = "none"),
#'     nrow = 1)
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

  idge <- .interaction_fdge(x, y, treat_lfc = treat_lfc, rerun = rerun, ...)
  if (is.null(idge[["result"]])) {
    # this happens if idge is null, too
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

  if (!is.null(idge[["result"]])) {
    ires <- tidy(idge[["result"]]) %>%
      select(feature_type, feature_id, {{stat.cols}})
    xystats <- left_join(xystats, ires, by = c("feature_type", "feature_id"))
    # put stats for interaction test up front, followed by *.x, *.y
    # xystats <- select(xystats, !!c(meta.cols, stat.cols), everything())
    xystats <- select(xystats, {{meta.cols}}, {{stat.cols}}, everything())
    out[["result"]] <- idge[["result"]]
  } else {
    # out[["result"]] <- NULL
  }

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

#' Decorates statistics table w/ individual/joint significance calls
#'
#' @noRd
#' @export
tidy.FacileTtestComparisonAnalysisResult <- function(
    x, max_padj_x = 0.1, min_logFC_x = NULL,
    max_padj_y = max_padj_x, min_logFC_y = min_logFC_x,
    labels = NULL, ...) {

  out <- x[["xystats"]]
  xres <- param(x, "x")
  yres <- param(x, "y")

  assert_number(max_padj_x, lower = 1e-6, upper = 1)
  assert_number(max_padj_y, lower = 1e-6, upper = 1)

  if (is.null(min_logFC_x)) {
    min_logFC_x <- param(xres, "treat_lfc")
    if (is.null(min_logFC_x)) min_logFC_x <- 1
  }
  assert_number(min_logFC_x, lower = 0, upper = Inf)
  if (is.null(min_logFC_y)) {
    min_logFC_y <- param(yres, "treat_lfc")
    if (is.null(min_logFC_y)) min_logFC_y <- 1
  }
  assert_number(min_logFC_y, lower = 0, upper = Inf)

  default.labels <- c(x = "x", y = "y", both = "both", none = "none")
  assert_character(labels, null.ok = TRUE, names = "unique")
  labels <- c(labels, default.labels)
  labels <- labels[!duplicated(names(labels))]
  assert_subset(c("x", "y", "both", "none"), names(labels))

  out <- mutate(
    out,
    sigclass = case_when(
      .data$padj.x <= .env$max_padj_x & .data$padj.y <= .env$max_padj_y ~ labels["both"],
      .data$padj.x <= .env$max_padj_x & .data$padj.y > .env$max_padj_y ~ labels["x"],
      .data$padj.y <= .env$max_padj_x & .data$padj.x > .env$max_padj_y ~ labels["y"],
      TRUE ~ labels["none"]))
  attr(out, "labels") <- labels
  out
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
                                                    features = NULL,
                                                    highlight = NULL,
                                                    facet = TRUE,
                                                    static = FALSE,
                                                    ...) {
  if (static) {
    ret <- sviz.FacileTtestComparisonAnalysisResult(x, max_padj = max_padj,
                                                    features = features,
                                                    highlight = highlight,
                                                    facet = facet, ...)
    return(ret)
  }

  hover <- c(
    # feature metadata
    "symbol", "feature_id", "meta",
    # padj from individual fdge tests
    "padj.x", "padj.y",
    # stats from interaction model, if it was run
    "logFC", "padj")
  d <- tidy(x)

  features <- extract_feature_id(features)
  highlight <- extract_feature_id(highlight)
  ids <- NULL

  if (!is.null(x[["dge"]])) {
    ids <- filter(d, padj.x <= max_padj | padj.y <= max_padj | padj <= max_padj)
    ids <- ids[["feature_id"]]
  } else {
    ids <- filter(d, padj.x <= max_padj | padj.y <= max_padj)[["feature_id"]]
  }

  dat <- filter(d, .data$feature_id %in% c(ids, features, highlight))

  if (length(highlight)) {
    dat[["highlight."]] <- ifelse(dat[["feature_id"]] %in% highlight,
                                  "fg", "bg")
    color_aes <- "highlight."
  } else {
    color_aes <- NULL
  }

  hover <- setdiff(intersect(hover, colnames(d)), "highlight.")
  fscatterplot(dat, c("logFC.x", "logFC.y"), hover = hover, webgl = TRUE,
               color_aes = color_aes, showlegend = FALSE, ...)
}

#' @noRd
#' @export
sviz.FacileTtestComparisonAnalysisResult <- function(x, max_padj = 0.1,
                                                     features = NULL,
                                                     highlight = NULL,
                                                     cor.method = "spearman",
                                                     cor.use = "complete.obs",
                                                     title = "DGE Comparison",
                                                     subtitle = NULL,
                                                     with_cor = TRUE, ...) {
  xdat <- tidy(x, max_padj_x = max_padj, max_pady_y = max_padj, ...)
  if (with_cor) {
    cors.all <- lapply(c("all", unique(xdat[["sigclass"]])), function(wut) {
      xs <- if (wut == "all") xdat else filter(xdat, .data$sigclass == .env$wut)
      xs <- filter(xs, !is.na(logFC.x) & !is.na(logFC.y))
      if (nrow(xs) == 0) return(NULL)
      ct <- suppressWarnings(
        cor.test(xs$logFC.x, xs$logFC.y, method = cor.method)
      )
      mutate(tidy(ct), sigclass = wut, n = nrow(xs))
    })
    cors <- mutate(bind_rows(cors.all),
                   label = sprintf("cor: %0.2f\nN: %d", estimate, n))
  }

  labels <- attr(xdat, "labels")[c("none", "both", "x", "y")]
  # labels <- factor(labels, labels)
  xdat[["sigclass"]] <- factor(xdat[["sigclass"]], unname(labels))

  cols.comp <- setNames(
    c("lightgrey", "darkgrey", "cornflowerblue", "orange"),
    labels)

  lims.square <- range(c(xdat$logFC.x, xdat$logFC.y))
  lims.square <- c(-1, 1) * (max(abs(lims.square)) + 0.1)

  if (is.null(subtitle) && !(labels["x"] == "xsig" || labels["y"] == "ysig")) {
    subtitle <- sprintf("%s vs %s", labels["x"], labels["y"])
  }

  lnone <- labels["none"]
  lboth <- labels["both"]

  gg.base <- xdat %>%
    ggplot2::ggplot(ggplot2::aes(x = logFC.x, y = logFC.y)) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = sigclass),
                        data = filter(xdat, .data$sigclass == .env$lnone)) +
    ggplot2::geom_point(ggplot2::aes(color = sigclass),
                        data = filter(xdat, .data$sigclass == .env$lboth)) +
    ggplot2::geom_point(ggplot2::aes(color = sigclass),
                        data = filter(xdat, !.data$sigclass %in% c(lnone, lboth))) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
    ggplot2::scale_color_manual(values = cols.comp) +
    ggplot2::labs(
      x = sprintf("log2FC %s", labels["x"]),
      y = sprintf("log2FC %s", labels["y"]),
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::xlim(lims.square) +
    ggplot2::ylim(lims.square)

  if (with_cor) {
    gg.main <- gg.base +
      ggplot2::geom_text(
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        data = filter(cors, sigclass == "all"))
  } else {
    gg.main <- gg.base
  }

  gg.facets <- gg.base +
    ggplot2::facet_wrap(~ sigclass) +
    ggplot2::ylab(NULL) +
    ggplot2::labs(title = NULL, subtitle = NULL)

  if (with_cor) {
    fcors <- filter(cors, sigclass != "all")
    fcors[["sigclass"]] <- factor(fcors[["sigclass"]],
                                  levels(xdat[["sigclass"]]))
    gg.facets <- gg.facets +
      ggplot2::geom_text(
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        data = filter(fcors))
  }

  out <- list(
    plot = gg.main,
    plot_facets = gg.facets,
    input_data = xdat,
    params = list())

  class(out) <- c("FacileTtestComparisonViz", "FacileStaticViz", "FacileViz")
  out
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
          "https://github.com/facilebio/FacileAnalysis/issues/19",
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

  args <- list(...)
  imodel <- flm_def(samples., icovariate, batch = ibatch,
                    block = args[["block"]], contrast. = contrast.)
  genes. <- unique(c(xres[["feature_id"]], yres[["feature_id"]]))

  ires <- tryCatch(
    fdge(imodel, features = genes.,
         method = param(x, "method"),
         assay_name = param(x, "assay_name"),
         with_sample_weights = param(x, "with_sample_weights"),
         treat_lfc = treat_lfc),
    error = function(x) NULL)

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
