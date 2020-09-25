#' @include fdge.R
NULL

#' @rdname fdge
#' @export
#'
#' @section Comparing DGE Results (interaction effect):
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
#' comp.hi <- tidy(dge.comp) %>%
#'   dplyr::group_by(interaction_group) %>%
#'   dplyr::slice(1:3) %>%
#'   dplyr::ungroup()
#' # Static visualization generates the main "4-way" plot, as well as the
#' # facets for each category.
#' sviz <- viz(dge.comp, labels = c(x = "CRC", y = "BLCA"),
#'             subtitle = "Tumor vs normal comparisons across indications",
#'             highlight = comp.hi)
#' # highlight some of them
#' s.hi <- sviz$input_data %>%
#'   dplyr::group_by(interaction_group) %>%
#'   dplyr::slice(1:3) %>%
#'   dplyr::ungroup()
#' if (requireNamespace("patchwork")) {
#'   patchwork::wrap_plots(
#'     sviz$plot + ggplot2::theme(legend.position = "bottom"),
#'     sviz$plot_facets + ggplot2::theme(legend.position = "none"),
#'     nrow = 1)
#'   viz(dge.comp, labels = c(x = "CRC", y = "BLCA"),
#'       color_quadrant = "darkgrey")$plot_facets
#' }
compare.FacileTtestAnalysisResult <- function(x, y,
                                              treat_lfc = param(x, "treat_lfc"),
                                              rerun = TRUE, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  # TODO: assert_comparable(x, y, ...)
  # TODO: assert_comparable.FacileDgeAnalysisResult <- function(...)
  assert_class(x, "FacileTtestAnalysisResult")
  assert_class(y, "FacileTtestAnalysisResult")
  fds. <- assert_class(fds(x), "FacileDataStore")

  # are there any joint features found in here?
  if (length(intersect(features(x)$feature_id, features(y)$feature_id)) == 0) {
    stop("No joint features found, comparison is not possible")
  }

  assert_flag(rerun)
  if (is.null(treat_lfc)) treat_lfc <- 0
  assert_number(treat_lfc, lower = 0)

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

  if (!isTRUE(idge[["with_stats"]])) {
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

  meta.cols <- c("feature_type", "feature_id", "symbol", "name", "meta")
  drop.cols <- c("seqnames", "start", "end", "strand", "effective_length",
                 "source")
  stat.cols <- setdiff(colnames(xres),  c(meta.cols, drop.cols))

  meta.cols <- intersect(meta.cols, jcols)
  xystats <- full_join(
    select(xres, {{meta.cols}}, {{stat.cols}}),
    select(yres, feature_type, feature_id, {{stat.cols}}),
    by = c("feature_type", "feature_id"))

  if (isTRUE(idge[["with_stats"]])) {
    ires <- tidy(idge[["result"]]) %>%
      select(feature_type, feature_id, {{stat.cols}})
  } else {
    ires <- idge[["result"]]
  }
  xystats <- left_join(xystats, ires, by = c("feature_type", "feature_id"))
  # put stats for interaction test up front, followed by *.x, *.y
  # xystats <- select(xystats, !!c(meta.cols, stat.cols), everything())
  take.cols <- intersect(c(meta.cols, stat.cols), colnames(xystats))
  xystats <- select(xystats, {{take.cols}}, everything())
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

#' @section Statistics Tables:
#' The stats table from the differential expression analysis can be retrieved
#' using the `tidy()` function. Depending upon the type of analysis run, the
#' exact columns of the returned table may differ.
#'
#' **Interaction Statistics**
#' Calling `tidy()` on an interaction test result
#' (`FacileTtestComparisonAnalysisResult`), returns the statistics for the
#' interaction test itself (if it was performed), as well as the statistics
#' for the individual DGE results that were run vai `compare(x, y)` to get the
#' interaction results itself. The columns of statistics related to the
#' individual tests will be suffixed with `*.x` and `*.y`, respectively.
#'
#' An `interaction_group` column will also be added to indicate what type of
#' statisticaly significance was found for each gene. The values in here can be:
#'
#' 1. `"both"`: this gene was statistically significant in both tests
#' 2. `"x"`: this gene was only significant in the `x` dge result
#' 3. `"y"`: this gene was only significant in the `y` dge result
#' 4. `"none"`: was not statistically significant in either test result
#'
#' The genes selected for significance from the input results `x` and `y` are
#' based on their `padj` and `logFC` values. The thresholds are tweaked by the
#' following parameters in the call to `tidy(compare(x,y))`:
#'
#' 1. `max_padj_(x|y)`: if not specified, defaults to `0.10`
#' 2. `min_logFC_(x|y)`: if not specified, we will take the value that was
#'    used in the `treat_lfc` parameters to `x` and `y` if those tests were
#'    run against a threshold, otherwise defaults to `1`.
#'
#' @rdname fdge
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

  .xsig <- out$padj.x <= max_padj_x & abs(out$logFC.x) >= min_logFC_x
  .ysig <- out$padj.y <= max_padj_y & abs(out$logFC.y) >= min_logFC_y

  out <- mutate(
    out,
    interaction_group = case_when(
       .xsig &  .ysig         ~ labels["both"],
       .xsig & !.ysig         ~ labels["x"],
      !.xsig &  .ysig         ~ labels["y"],
      TRUE                    ~ labels["none"])) %>%
    select(interaction_group, everything())
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

#' If `interacive` is `TRUE` (default), then this visualizaiton will drop the
#' interactive points by setting `insignificant = "drop"`. Static plots keep the
#' insignificant points. If you want to change this behavior, set `interactive`
#' and `insignificant` as you please.
#' @noRd
#' @export
viz.FacileTtestComparisonAnalysisResult <- function(
    x, max_padj = 0.1, features = NULL, highlight = NULL,
    color_quadrant = NULL, color_highlight = "red",
    cor.method = "spearman", title = "DGE Comparison",
    subtitle = NULL, with_cor = TRUE, interactive = TRUE,
    insignificant = if (interactive) "drop" else "points",
    facets_nrow = 2, ...) {

  xdat <- tidy(x, max_padj_x = max_padj, max_pady_y = max_padj, ...)
  labels <- attr(xdat, "labels")[c("none", "both", "x", "y")]
  insignificant <- match.arg(insignificant, c("drop", "points"))

  default.quadrant.cols <- setNames(
    c("lightgrey", "darkgrey", "cornflowerblue", "orange"),
    labels)
  if (is.null(color_quadrant) || length(color_quadrant) == 0L) {
    color_quadrant <- default.quadrant.cols
  } else {
    assert_character(color_quadrant, min.len = 1)
    if (length(color_quadrant) == 1L) {
      color_quadrant <- rep(color_quadrant, length(default.quadrant.cols))
      names(color_quadrant) <- names(default.quadrant.cols)
    }
    color_quadrant <- ifelse(
      is.na(color_quadrant[names(default.quadrant.cols)]),
      default.quadrant.cols,
      color_quadrant)
  }

  if (insignificant == "drop") {
    xdat <- filter(xdat, interaction_group != labels["none"])
    labels <- labels[-1]
    color_quadrant <- color_quadrant[-1]
  }

  xdat[["interaction_group"]] <- factor(xdat[["interaction_group"]],
                                        unname(labels))

  if (with_cor) {
    cor.quadrants <- c("all", levels(xdat[["interaction_group"]]))
    cors.all <- lapply(cor.quadrants, function(wut) {
      if (wut == "all") {
        xs <- xdat
      } else {
        xs <- filter(xdat, .data$interaction_group == .env$wut)
      }
      xs <- filter(xs, !is.na(logFC.x) & !is.na(logFC.y))
      if (nrow(xs) == 0) return(NULL)
      ct <- suppressWarnings(
        cor.test(xs$logFC.x, xs$logFC.y, method = cor.method)
      )
      mutate(tidy(ct), interaction_group = wut, n = nrow(xs))
    })
    cors <- mutate(bind_rows(cors.all),
                   label = sprintf("cor: %0.2f\nN: %d", estimate, n))
  }

  lims.square <- range(c(xdat$logFC.x, xdat$logFC.y))
  lims.square <- c(-1, 1) * (max(abs(lims.square)) + 0.1)

  if (is.null(subtitle) && !(labels["x"] == "xsig" || labels["y"] == "ysig")) {
    subtitle <- sprintf("%s vs %s", labels["x"], labels["y"])
  }

  lnone <- labels["none"]
  lboth <- labels["both"]

  xdat$padj.min <- pmin(xdat$padj, xdat$padj.x, xdat$padj.y)
  xdat <- arrange(xdat, interaction_group, desc(padj.min))

  gg.base <- xdat %>%
    ggplot2::ggplot(ggplot2::aes(x = logFC.x, y = logFC.y)) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed")
  # if (insignificant == "density") {
  #   gg.base <- gg.base +
  #     ggplot2::geom_density2d(
  #       alpha = 0.5,
  #       data = filter(xdat, anti_join(xdat, xdat.points, by = "feature_id")))
  # }
  gg.base <- gg.base +
    suppressWarnings({
      ggplot2::geom_point(
        ggplot2::aes(color = interaction_group, text = symbol))
    })

  gg.base <- gg.base +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dotted") +
    ggplot2::scale_color_manual(values = color_quadrant) +
    ggplot2::labs(
      x = sprintf("log2FC %s", labels["x"]),
      y = sprintf("log2FC %s", labels["y"]),
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::xlim(lims.square) +
    ggplot2::ylim(lims.square)

  if (!is.null(highlight)) {
    fids <- extract_feature_id(highlight)
    highlight <- filter(xdat, feature_id %in% .env$fids)
    if (nrow(highlight)) {
      gg.base <- gg.base +
        suppressWarnings({
          ggplot2::geom_point(
            ggplot2::aes(text = symbol),
            data = highlight,
            color = color_highlight)
        })
    }
  }

  if (with_cor) {
    gg.main <- gg.base +
      ggplot2::geom_text(
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        data = filter(cors, interaction_group == "all"))
  } else {
    gg.main <- gg.base
  }

  gg.facets <- gg.base +
    ggplot2::facet_wrap(~ interaction_group, nrow = facets_nrow) +
    ggplot2::ylab(NULL) +
    ggplot2::labs(title = NULL, subtitle = NULL)

  if (with_cor) {
    fcors <- filter(cors, interaction_group != "all")
    fcors[["interaction_group"]] <- factor(fcors[["interaction_group"]],
                                           levels(xdat[["interaction_group"]]))
    gg.facets <- gg.facets +
      ggplot2::geom_text(
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        data = fcors)
  }

  out <- list(
    plot = gg.main,
    plot_facets = gg.facets,
    input_data = xdat,
    correlation = cors,
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
.interaction_fdge <- function(x, y, treat_lfc = NULL, rerun = FALSE, ...,
                              .run_interaction = TRUE) {
  # If these results aren't from the same FacileDataStore, get outta here
  # NOTE: This is not really a robust way to compare if two fds are the same
  xmod <- model(x)
  ymod <- model(y)
  xres <- tidy(x)
  yres <- tidy(y)

  covariate <- param(xmod, "covariate")

  concordant.analysis <- covariate == param(ymod, "covariate") &&
    param(x, "assay_name") == param(y, "assay_name") &&
    param(x, "method") == param(y, "method")

  # Calculate estimated logFC's between the *.x and *.y logFC in case we don't
  # make it all the way through the formal interaction analysis
  ires.tmp <- xres %>%
    full_join(yres, by = c("feature_id", "feature_type", "name")) %>%
    transmute(feature_id, feature_type, name, logFC = logFC.x - logFC.y,
              pval = NA_real_, padj = NA_real_)

  out <- list(result = ires.tmp, x = x, y = y, rerun = FALSE,
              with_stats = FALSE)
  if (!concordant.analysis) {
    warning("Analyses are not concordant. Interaction stats will not be ",
            "generated, but delta logFC's will be provided")
    return(out)
  }
  if (!.run_interaction) {
    return(out)
  }

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

  if (!is.null(ires)) {
    out[["result"]] <- ires
    out[["with_stats"]] <- TRUE
    rerun <- rerun && !setequal(xres[["feature_id"]], genes.)

    if (rerun) {
      x <- fdge(xmod, features = genes., method = param(x, "method"),
                assay_name = param(x, "assay_name"),
                with_sample_weights = param(x, "with_sample_weights"))
      y <- fdge(ymod, features = genes., method = param(y, "method"),
                assay_name = param(y, "assay_name"),
                with_sample_weights = param(y, "with_sample_weights"))
      out[["rerun"]] <- TRUE
      out[["x"]] <- x
      out[["y"]] <- y
    }
  }

  out
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
