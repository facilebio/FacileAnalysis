#' Construct a bioconductor classed object from an analysis.
#'
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
                                                method = NULL, features = NULL,
                                                filter = "default",
                                                filter_universe = NULL,
                                                filter_require = NULL,
                                                with_sample_weights = FALSE,
                                                weights = NULL, block = NULL,
                                                prior_count = 0.1, ...) {
  assert_class(x, "FacileLinearModelDefinition")
  si <- assert_class(x$covariates, "facile_frame")
  .fds <- assert_class(fds(x), "FacileDataStore")
  if (is.null(assay_name)) assay_name <- default_assay(.fds)

  # Because this returns a Bioconductor object, we will store the expected
  # (internal) facile bits in the `ifacile` attribute
  ifacile <-  list(
    params = list(
      assay_name = assay_name,
      method = method,
      features = features,
      filter = filter,
      with_sample_weights = with_sample_weights,
      prior_count = prior_count))

  messages <- character()
  warnings <- character()
  errors <- character()

  out <- try({
    stop('biocbox did not materialize, see `attr(this, "ifacile")$errors`')
  }, silent = TRUE)

  on.exit({
    ifacile[["messages"]] <- messages
    ifacile[["warnings"]] <- warnings
    ifacile[["errors"]] <- errors
    attr(out, "facile") <- ifacile
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

  if (isTRUE(filter)) filter <- "default"
  if (is.null(filter)) filter <- FALSE # originally NULL was equal to FALSE
  if (!(test_string(filter) || isFALSE(filter))) {
    errors <- c(
      glue("Invalid type for `filter` argument (`{class(filter)}`).",
           "Only a string or logical flag is allowed"))
    return(out)
  }

  if (is.null(method)) {
    method <- valid_methods[["dge_method"]][1L]
  }
  if (!method %in% valid_methods[["dge_method"]]) {
    default_method <- valid_methods[["dge_method"]][1L]
    msg <- glue("Requested dge_method `{method}` not found, using ",
                "`{default_method}` instead")
    warnings <- c(warnings, msg)
    method <- default_method
  }

  ifacile[["params"]][["method"]] <- method

  method.info <- filter(valid_methods, dge_method == method)
  bb.all <- biocbox(si, class = method.info[["bioc_class"]],
                    assay_name = assay_name,
                    features = features, sample_covariates = FALSE)
  if (is.null(features)) {
    # Only execute a filtering strategy if the caller did not explicitly request
    # features.
    bb <- .filter_features(bb.all, x, ainfo, method.info$default_filter,
                           filter = filter, filter_universe = filter_universe,
                           filter_require = filter_require, ...)
  } else {
    bb <- bb.all
    if (is(bb, "DGEList") && all(bb[["samples"]][["norm.factors"]] == 1)) {
      bb <- edgeR::calcNormFactors(bb)
    }
  }

  des <- design(x)
  if (!setequal(rownames(des), colnames(bb))) {
    # This can be triggered when the samples the design was defined on do not
    # have molecular data from the `assay_name` we are testing against.
    warning("Reducing size of samples", immediate. = TRUE)
    # errors <- c(errors, "rownames of design(x) does not match biocbox")
    # return(out)
  }
  if (!isTRUE(all.equal(rownames(des), colnames(bb)))) {
    bb <- bb[,rownames(des)]
  }

  # assert rownames(des) == colnames(bb)
  if (method == "edgeR-qlf") {
    out <- suppressWarnings(edgeR::estimateDisp(bb, des, robust = TRUE))
  } else if (method == "voom") {
    # out <- .voomLmFit(bb, des, block = param(x, "block"),
    #                   prior.weights = weights,
    #                   sample.weights = with_sample_weights)

    if (with_sample_weights) {
      out <- .voomWithQualityWeights(bb, des, block = param(x, "block"), ...)
    } else {
      out <- .voom(bb, des, block = param(x, "block"), ...)
    }

    if (!is.null(weights)) {
      warning("Can't estimate sample weights and include prior weights. ",
              "The prior `weights` you provided have been set to NULL",
              immediate. = TRUE)
      weights <- NULL
    }
  } else if (method %in% c("limma-trend", "limma")) {
    if (is(bb, "DGEList")) {
      prior_count <- max(0.1, prior_count)
      out <- new(
        "EList",
        list(E = edgeR::cpm(bb, log = TRUE, prior.count = prior_count),
             genes = bb[["genes"]], targets = bb[["samples"]]))
    } else {
      out <- bb
    }
    stopifnot(is(out, "EList"))
    out[["design"]] <- des
    if (with_sample_weights) {
      out[["weights"]] <- limma::arrayWeights(out, out[["design"]])
    }
  } else {
    stop("How did we get here?")
  }

  if (!is.null(weights)) {
    out <- .add_observation_weights(out, weights, ...)
  }

  out
}

#' @noRd
.add_observation_weights <- function(x, weights, ...) {
  if (!is(x, "EList")) {
    warning("observation weights only supported for ELists. They are being ",
            "ignored.", immediate. = TRUE)
    return(x)
  }

  if (!is.null(weights)) {
    if (!is.null(x[["weights"]])) {
      warning("weights already calculated, but are being replaced",
              immediate. = TRUE)
    }
    assert_multi_class(weights, c("data.frame", "tibble"))
    assert_subset(c("dataset", "sample_id",  "feature_id", "weight"),
                  colnames(weights))

    weights <- mutate(weights, skey = paste(dataset, sample_id, sep = "__"))

    ww <- select(weights, skey, feature_id, weight)
    ww <- pivot_wider(ww, names_from = skey, values_from = weight)
    weights <- as.matrix(select(ww, -feature_id))
    rownames(weights) <- ww[["feature_id"]]
    weights <- weights[rownames(x), colnames(x)]

    na.w <- is.na(weights)
    if (any(na.w)) {
      mean.w <- mean(weights[!na.w])
      msg <- sprintf("%0.2f NA values in weights, replacing with mean %0.2f",
                     sum(na.w) / length(na.w), mean.w)
      warning(msg)
      weights[na.w] <- mean.w
    }
    x[["weights"]] <- weights
  }

  x
}

# Feature abundance filtering helper functions ---------------------------------

#' @noRd
#' @param default.filter the default feature filtering option for the dge method
#'   and assay_type we are prepping for. (TRUE/FALSE)
.filter_features <- function(x, model, ainfo, default.filter, filter,
                             filter_universe, filter_require,
                             update_normfactors = NULL,
                             ...) {
  assert_flag(default.filter)
  if (is(x, "DGEList")) {
    if (is.null(update_normfactors)) {
      update_normfactors <- all(x[["samples"]][["norm.factors"]] == 1)
    }
    if (isTRUE(update_normfactors)) {
      x <- edgeR::calcNormFactors(x)
    }
  }

  filter_require <- extract_feature_id(filter_require)
  filter_universe <- extract_feature_id(filter_universe)
  if (!is.null(filter_universe)) {
    filter_universe <- unique(c(filter_universe, filter_require))
  }

  do.filterByExpr <- !isFALSE(filter) &&
    (isTRUE(filter) || (isTRUE(filter == "default") && isTRUE(default.filter)))

  if (do.filterByExpr) {
    if (!ainfo[["assay_type"]] %in% c("rnaseq", "umi", "pseudobulk")) {
      warning("Automatic filtering only implemented for `rnaseq` and `umi` ",
              "assay types for now.\nFiltering step skipped",
              immediate. = TRUE)
      return(x)
    }
    des.matrix <- design(model)[colnames(x),,drop=FALSE]
    # keep only the columns in the design matrix used for filtering that
    # correspond to the main groups. If this is a ttest, its
    # des.matrix[, design$text_covs]. If anova, then we also include the
    # intercept column
    filter.cols <- model$test_covs
    if (is_anova(model)) {
      filter.cols <- c(
        grep("(Intercept)", colnames(des.matrix), ignore.case = TRUE),
        filter.cols)
    }
    x <- .filter_count_assay(x, des.matrix, filter.cols,
                             filter_universe, filter_require, ...)
  }
  x
}

#' @noRd
.filter_count_assay <- function(x, des.matrix, filter_design_columns,
                                filter_universe, filter_require,
                                # default params for edgeR::filterByExpr
                                filter_min_count = 10,
                                filter_min_total_count = 15,
                                ...) {
  if (!is(x, "DGEList")) {
    warning("Automatic filtering only implemented for DGELists for now",
            immediate. = TRUE)
    return(x)
  }
  dmatrix <- des.matrix[, filter_design_columns, drop = FALSE]
  dmatrix <- dmatrix[rowSums(dmatrix) > 0,,drop = FALSE]
  if (is.character(filter_universe)) {
    x <- x[rownames(x) %in% filter_universe,]
  }
  keep <- edgeR::filterByExpr(x[, rownames(dmatrix)], dmatrix,
                              min.count = filter_min_count,
                              min.total.count = filter_min_total_count, ...)
  if (is.character(filter_require)) {
    keep <- keep | rownames(x) %in% filter_require
  }
  keep_fraction <- mean(keep)
  keep_n <- sum(keep)

  if (keep_fraction < 0.50) {
    msg <- glue("Only {format(keep_fraction * 100, digits = 4)}% ",
                "({keep_n}) of features are retained after filtering.")
    # warnings <- c(warnings, msg)
    warning(msg, immediate. = TRUE)
  }

  x <- x[keep,,keep.lib.sizes = FALSE]
  suppressWarnings(edgeR::calcNormFactors(x)) # partial match of `p` to `probs`
}

