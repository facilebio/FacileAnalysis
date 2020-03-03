# This was the original biocbox creation and filtering stuff before we pushed
# the biocbox S3 generic up to FacileData and now have
# biocbox.FacileLinearModelDefinition() return a "naked" bioconductor container

#' Something like a factor function for a "biocbox"
#'
#' All the argument checking and whatever else has been done upstream, these
#' helper functions are just executing under the assumption that everything
#' is kosher.
#'
#' We'll come back to make this pretty later (#!refactor)
#'
#' @noRd
.biocbox_create <- function(xsamples, assay_name, assay_type,
                            design, filter, method, with_sample_weights,
                            prior_count, ...) {
  assert_class(design, "FacileLinearModelDefinition")

  amethods <- fdge_methods(assay_type)
  if (nrow(amethods) == 0) {
    stop(paste("DGE analysis not implemented for this assay_type: ",
               assay_type))
  }
  xmethod <- filter(amethods, dge_method == method)
  if (nrow(xmethod) != 1L) {
    stop("Unexpected state, method parameter not legal: ", method)
  }
  bioc.class <- xmethod[["bioc_class"]]
  if (bioc.class == "DGEList") {
    create <- .biocbox_create_DGEList
  } else if (bioc.class == "EList") {
    create <- .biocbox_create_EList
  } else {
    stop("Unexpxected state, ildefined bioc.class for method: ", bioc.class)
  }

  bbox <- create(xsamples, assay_name = assay_name,
                 assay_type = assay_type,
                 design = design, filter = filter,
                 method = method,
                 with_sample_weights = with_sample_weights,
                 prior_count = prior_count, ...)
  bbox
}

#' @noRd
.get_DGEList <- function(xsamples, design, assay_name, assay_type, method,
                         filter, filter_universe, filter_require, ...) {
  if (test_multi_class(filter_require, c("data.frame", "tbl"))) {
    filter_require <- filter_require[["feature_id"]]
  }
  if (!is.null(filter_universe)) {
    if (test_multi_class(filter_universe, c("data.frame", "tbl"))) {
      filter_universe <- filter_universe[["feature_id"]]
      if (is.character(filter_require)) {
        filter_universe <- unique(c(filter_universe, filter_require))
      }
    } else {
      warning("Unknown parameter type for filter_universe, ignoring ...")
      filter_universe <- NULL
    }
  }

  amethod <- filter(fdge_methods(assay_type), dge_method == method)

  do.filterByExpr <- nrow(amethod) == 1L &&
    test_string(filter) &&
    filter == "default" &&
    isTRUE(amethod[["default_filter"]])

  if (!do.filterByExpr) {
    if (test_multi_class(filter, c("data.frame", "tbl"))) {
      filter_universe <- filter[["feature_id"]]
    } else if (test_string(filter) && filter == "default") {
      filter_universe <- NULL
    } else if (test_character(filter)) {
      filter_universe <- filter
    }
    if (!is.character(filter_universe) && !is.null(filter_universe)) {
      stop("filter argument is of illegal type: ", class(filter)[1L])
    }
    filter.cols <- NULL
  }

  update_libsizes <- FALSE
  update_normfactors <- FALSE

  if (is.character(filter_universe)) {
    filter_universe <- unique(c(filter_universe, filter_require))
  }

  y <- as.DGEList(xsamples, assay_name = assay_name,
                  covariates = xsamples,
                  feature_ids = filter_universe,
                  update_libsizes = update_libsizes,
                  update_normfactors = update_normfactors)

  if (do.filterByExpr) {
    des.matrix <- design(design)[colnames(y),,drop=FALSE]
    # keep only the columns in the design matrix used for filtering that
    # correspond to the main groups. If this is a ttest, its
    # des.matrix[, design$text_covs]. If anova, then we also include the
    # intercept column
    filter.cols <- design$test_covs
    if (is.anova(design)) {
      filter.cols <- c(
        grep("(Intercept)", colnames(des.matrix), ignore.case = TRUE),
        filter.cols)
    }
  }

  list(
    y = y,
    filter_run = do.filterByExpr,
    filter_require = filter_require,
    filter_design_columns = filter.cols)
}

#' The behavior of the `filter*` parameters are currently described in the
#' [fdge()] help page.
#'
#' @noRd
#' @importFrom edgeR filterByExpr calcNormFactors estimateDisp
#' @importFrom limma arrayWeights voom voomWithQualityWeights
.biocbox_create_DGEList <- function(xsamples, assay_name, assay_type,
                                    design, filter, method, with_sample_weights,
                                    prior_count = 2,
                                    # default params for edgeR::filterByExpr
                                    filter_min_count = 10,
                                    filter_min_total_count = 15,
                                    # Additional filter params,
                                    filter_universe = NULL,
                                    filter_require = NULL,
                                    # minimum gene count to reset lib.size and
                                    # norm factors
                                    min_feature_count_tmm = 500, ...) {
  assert_class(design, "FacileLinearModelDefinition")
  if (is.null(prior_count)) prior_count <- 2
  warnings <- character()

  dat <- .get_DGEList(xsamples, design, assay_name = assay_name,
                      assay_type = assay_type, method = method,
                      filter = filter,
                      filter_universe = filter_universe,
                      filter_require = filter_require)
  y <- dat[["y"]]
  des.matrix <- design(design)[colnames(y),,drop=FALSE]

  if (dat[["filter_run"]]) {
    dmatrix <- des.matrix[, dat[["filter_design_columns"]], drop = FALSE]
    dmatrix <- dmatrix[rowSums(dmatrix) > 0,,drop = FALSE]
    keep <- filterByExpr(y[, rownames(dmatrix)], dmatrix,
                         min.count = filter_min_count,
                         min.total.count = filter_min_total_count, ...)
    if (is.character(dat[["filter_require"]])) {
      keep <- keep | rownames(y) %in% dat[["filter_require"]]
    }
    keep_fraction <- mean(keep)
    keep_n <- sum(keep)

    if (keep_fraction < 0.50 * nrow(y)) {
      msg <- glue("Only {format(keep_fraction * 100, digits = 4)}% ",
                  "({keep_n}) of features are retained after filtering.")
      warnings <- c(warnings, msg)
    }

    y <- y[keep,,keep.lib.sizes = FALSE]
    y <- suppressWarnings(calcNormFactors(y)) # partial match of `p` to `probs`
  } else {
    # update lib.size and normfactors if enough genes here
    if (nrow(y) > min_feature_count_tmm) {
      y <- y[rep(TRUE, nrow(y)),,keep.lib.sizes = FALSE]
      y <- suppressWarnings(calcNormFactors(y)) # partial match of `p` to `probs`
    }
  }

  y$design <- des.matrix

  if (method == "edgeR-qlf") {
    out <- suppressWarnings(estimateDisp(y, y[["design"]], robust = TRUE))
  } else if (method == "voom") {
    if (with_sample_weights) {
      out <- voomWithQualityWeights(y, y[["design"]], ...)
    } else {
      out <- .voom_dots(y, y$design, ...)
    }
  } else if (method %in% c("limma-trend", "limma")) {
    elist <- list()
    elist[["E"]] <- edgeR::cpm(y, log = TRUE, prior.count = prior_count)
    elist[["genes"]] <- y[["genes"]]
    elist[["targets"]] <- y[["samples"]]
    elist <- new("EList", elist)
    elist[["design"]] <- y[["design"]]
    if (with_sample_weights) {
      elist <- arrayWeights(elist, elist[["design"]])
    }
    out <- elist
  } else {
    stop("How did we get here?")
  }

  dropped <- setdiff(rownames(dat[["y"]]), rownames(out))
  attr(out, "warnings") <- warnings
  out
}

#' @noRd
#' @importFrom stats hat
#' @importFrom limma arrayWeights
.biocbox_create_EList <- function(xsamples, assay_name, assay_type,
                                  design, filter, method, with_sample_weights,
                                  prior_count = 0.2,
                                  filter_min_expr = 1,
                                  filter_universe = NULL,
                                  filter_require = character(), ...) {
  if (is.null(prior_count)) prior_count <- 0.25

  # assay_type is one of c("normcounts","lognorm")
  # Note that as.DGEList always assembles the assay matrix with
  # normalized = FALSE

  dat <- .get_DGEList(xsamples, design, assay_name = assay_name,
                      assay_type = assay_type, method = method,
                      filter = filter,
                      filter_universe = filter_universe,
                      filter_require = filter_require)
  y <- dat[["y"]]
  des.matrix <- design(design)[colnames(y),,drop=FALSE]

  e <- y[["counts"]]
  if (assay_type %in% c("normcounts")) {
    e <- log2(e + prior_count)
  }

  # Remove genes according to `filter` specificaiton
  if (dat[["filter_run"]]) {
    # if (assay_type == "lognorm") {
    #   min.expr <- log2(2)
    # } else {
    #   min.expr <- 1
    # }
    min.expr <- filter_min_expr
    dmatrix <- des.matrix[, dat[["filter_design_columns"]], drop = FALSE]
    dmatrix <- dmatrix[rowSums(dmatrix) > 0,,drop = FALSE]

    # Code below taken from edgeR:::filterByExpr.default function
    h <- hat(dmatrix)
    min.samples <- 1 / max(h)
    if (min.samples > 10) {
      min.samples <- 10 + (min.samples - 10) * 0.7
    }
    tol <- 1e-14
    keep <- rowSums(e[, rownames(dmatrix),drop = FALSE] >= min.expr) >= (min.samples - tol)
    # end filterByExpr block

    if (is.character(dat[["filter_require"]])) {
      keep <- keep | rownames(y) %in% dat[["filter_require"]]
    }
    keep_fraction <- mean(keep)
    keep_n <- sum(keep)

    if (keep_fraction < 0.50 * nrow(y)) {
      msg <- glue("Only {format(keep_fraction * 100, digits = 4)}% ",
                  "({keep_n}) of features are retained after filtering.")
      warnings <- c(warnings, msg)
    }
  } else {
    keep <- rep(TRUE, nrow(e))
  }

  elist <- list()
  elist[["E"]] <- e[keep,,drop = FALSE]
  elist[["genes"]] <- y[["genes"]][keep,,drop = FALSE]
  elist[["targets"]] <- y[["samples"]]
  elist <- new("EList", elist)
  elist[["design"]] <- des.matrix

  if (with_sample_weights) {
    elist <- arrayWeights(elist, elist[["design"]])
  }

  dropped <- setdiff(rownames(dat[["y"]]), rownames(elist))

  elist
}
