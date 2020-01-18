#' Materializes the bioc assay container to use to run a dge test.
#'
#' Bioconductor assay containers need to be materialized for certain analyses,
#' like differential gene expression. They also may need to be materialized
#' from a result of something generated within some FacileAnalysis.
#'
#' @export
#' @param x The FacileAnalysisResult
biocbox <- function(x, ...) {
  UseMethod("biocbox", x)
}

#' @noRd
#' @export
result.BiocBox <- function(x, ...) {
  x[["biocbox"]]
}

#' @noRd
#' @export
design.BiocBox <- function(x, ...) {
  box <- result(x)
  if (is.null(box)) {
    stop("The bioc assay container is not found in the BiocBox")
  }
  if (is(box, "DGEList") || is(box, "EList")) {
    out <- box[["design"]]
  } else {
    stop("design.BiocBox not implemented for '", class(box)[1L], "' containers")
  }

  if (!is.matrix(out)) {
    stop("Error extracting design matrix")
  }

  out
}

#' @noRd
#' @export
features.BiocBox <- function(x, ...)  {
  y <- result(x)
  if (is.null(y)) {
    stop("The bioc assay container is not found in the BiocBox")
  }
  if (test_multi_class(y, c("DGEList", "EList"))) {
    out <- y[["genes"]]
  } else {
    stop(class(y)[1L], " not yet handled")
  }

  if (is.null(out[["name"]]) && !is.null(out[["symbol"]])) {
    out[["name"]] <- out[["symbol"]]
  }

  as.tbl(out)
}

#' @section Linear Model Definitions:
#' This function accepts a model defined using using [flm_def()] and
#' creates the appropriate Bioconductor assay container to test the model
#' given the `assay_name` and dge `method` specified by the user.
#'
#' This function currently supports retrieving data and whipping it into
#' a DGEList (for count-like data) and an EList for data that can be analyzed
#' with one form limma or another.
#'
#' Assumptions on different `assay_type` values include:
#'
#' * `rnaseq`: assumed to be "vanilla" bulk rnaseq gene counts
#' * `umi`: data from bulk rnaseq, UMI data, like quantseq
#' * `tpm`: TPM values. These will be `log2(TPM + prior_count)` transformed,
#'          then differentially tested using the limma-trended pipeline
#'
#' TODO: support affymrna, affymirna, etc. assay types
#'
#' The "filter" parameters are described in the [fdge()] function for now.
#'
#' @export
#' @rdname biocbox
#'
#' @param sample_info a `facile_frame` that enumerates the samples to fetch
#'   data for, as well as the covariates used in downstream analysis
#' @param assay_name the name of the assay to pull data for
#' @param method the name of the dge method that will be used. This will dictate
#'   the post-processing of the data
#' @param filter A filtering policy to remove unintereesting genes.
#'   If `"default"` (which is the default), then [edgeR::filterByExpr()] is
#'   used if we are materializing a `DGEList`, otherwise lowly expressed
#'   features are removed in a similarly "naive" manner. This can,
#'   alternatively, be a character vector that holds the names of the features
#'   that should be kept. Default value: `"default"`.
#' @param with_sample_weights Some methods that leverage the limma pipeline,
#'   like `"voom"`, `"limma"`, and `"limma-trend"` can leverage sample (array)
#'   quality weights to downweight outlier samples. In the case of
#'   `method == "voom"`, we use [limma::voomWithQualityWeights()], while the
#'   rest use [limma::arrayWeights()]. The choice of `method` determines which
#'   sample weighting function to sue. Defaults to `FALSE`.
#' @param prior_count The pseudo-count to add to count data. Used primarily
#'   when running the `limma-trend` method on count (RNA-seq) data.
#' @param ... passed down to internal modeling and filtering functions.
#' @return a DGEList or EList with assay data in the correct place, and all of
#'   the covariates in the `$samples` or `$targerts` data.frame that are requied
#'   to test the model in `mdef`.
biocbox.FacileLinearModelDefinition <- function(x, assay_name = NULL,
                                                method = NULL,
                                                dge_methods = NULL,
                                                filter = "default",
                                                with_sample_weights = FALSE,
                                                prior_count = NULL, ...) {
  assert_class(x, "FacileLinearModelDefinition")
  si <- assert_class(x$covariates, "facile_frame")
  .fds <- assert_class(fds(x), "FacileDataStore")
  if (is.null(assay_name)) assay_name <- default_assay(.fds)

  out <- list(
    biocbox = NULL,
    params = list(
      assay_name = assay_name,
      method = method,
      filter = filter,
      with_sample_weights = with_sample_weights,
      prior_count = prior_count))
  class(out) <- "BiocBox"

  messages <- character()
  warnings <- character()
  errors <- character()

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    return(out)
  })

  ainfo <- assay_info(.fds, assay_name)
  assay_type <- ainfo[["assay_type"]]
  valid_methods <- fdge_methods(assay_type)
  if (nrow(valid_methods) == 0) {
    errors <- paste("DGE analysis not implemented for this assay_type: ",
                    assay_type)
    return(out)
  }

  if (is.factor(filter)) {
    filter <- as.character(filter)
  }
  if (!is.null(filter)) {
    if (test_multi_class(filter, c("data.frame", "tbl"))) {
      filter <- filter[["feature_id"]]
    }
    filter.kosher <- test_character(filter, min.len = 2) ||
      (test_string(filter) && filter == "default")
    if (!filter.kosher) {
      errors <- c(
        glue("Invalid `filter` value (`{filter}`). The only valid string value ",
             "the `filter` argument is 'default'"))
      return(out)
    }
  }

  if (is.null(dge_methods)) {
    dge_methods <- fdge_methods(ainfo[["assay_type"]])
  }
  if (is.null(method)) {
    method <- dge_methods[["dge_method"]][1L]
  }

  if (!method %in% dge_methods[["dge_method"]]) {
    default_method <- dge_methods[["dge_method"]][1L]
    msg <- glue("Requested dge_method `{method}` not found, using ",
                "`{default_method}` instead")
    warnings <- c(warnings, msg)
    method <- default_method
  }

  out[["biocbox"]] <- .biocbox_create(si, assay_name = assay_name,
                                      assay_type = assay_type,
                                      design = x, filter = filter,
                                      method = method,
                                      with_sample_weights = with_sample_weights,
                                      prior_count = prior_count, ...)
  out[["dge_method"]] <- method
  out
}

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
.get_DGEList <- function(xsamples, design, assay_name, assay_type,
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

  amethod <- fdge_methods(assay_type)

  do.filterByExpr <- test_string(filter) &&
    filter == "default" &&
    isTRUE(amethod[["default_filter"]])

  if (!do.filterByExpr) {
    if (test_multi_class(filter, c("data.frame", "tbl"))) {
      filter_universe <- filter[["feature_id"]]
    } else if (test_string(filter) && filter == "default") {
      filter_universe <- NULL
    }
    if (!is.character(filter_universe) && !is.null(filter_universe)) {
      stop("filter argument is of illegal type: ", class(filter)[1L])
    }
    filter.cols <- NULL
  }

  update_libsizes <- FALSE
  update_normfactors <- FALSE

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
                      assay_type = assay_type,
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

#' Defines a voom method that can accept dots (and toss them)
#'
#' @noRd
#' @importFrom limma voom
.voom_dots <- function(counts, design = NULL, lib.size = NULL,
                       normalize.method = "none",  block = NULL,
                       correlation = NULL, weights = NULL, span = 0.5,
                       plot = FALSE, save.plot = TRUE, ...) {
  voom(counts, design, lib.size = lib.size, normalize.method = normalize.method,
       block = block, correlation = correlation, weights = weights, span = span,
       plot = plot, save.plot = save.plot)
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
                      assay_type = assay_type,
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
