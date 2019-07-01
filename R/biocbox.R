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

#' @section Linear Model Definitions:
#' This function accepts a model defined using using [fdge_model_def()] and
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
biocbox.FacileDgeModelDefinition <- function(x, assay_name = NULL,
                                             method = NULL,
                                             dge_methods = NULL,
                                             filter = "default",
                                             with_sample_weights = FALSE,
                                             prior_count = NULL, ...) {
  assert_class(x, "FacileDgeModelDefinition")
  si <- assert_class(x$covariates, "facile_frame")
  .fds <- assert_class(fds(x), "FacileDataStore")
  if (is.null(assay_name)) assay_name <- default_assay(.fds)

  out <- list(biocbox = NULL)
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
  if (!assay_type %in% c("rnaseq", "umi", "normcounts", "lognorm")) {
    errors <- paste("DGE analysis not implemented for this assay_type: ",
                    assay_type)
    return(out)
  }

  if (test_string(filter) && filter != "default") {
    errors <- c(
      glue("Invalid `filter` value (`{filter}`). The only valid string value ",
            "the `filter` argument is 'default'"))
    return(out)
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
  assert_class(design, "FacileDgeModelDefinition")

  if (assay_type %in% c("rnaseq", "isoseq", "umi")) {
    create <- .biocbox_create_DGEList
  } else if (assay_type %in% c("normcounts", "lognorm")) {
    create <- .biocbox_create_EList
  } else {
    stop(paste("DGE analysis not implemented for this assay_type: ",
               assay_type))
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
#' @importFrom edgeR filterByExpr calcNormFactors estimateDisp
#' @importFrom limma arrayWeights voom voomWithQualityWeights
.biocbox_create_DGEList <- function(xsamples, assay_name, assay_type,
                                    design, filter, method, with_sample_weights,
                                    prior_count = 2,
                                    # default params for edgeR::filterByExpr
                                    filter_min_count = 10,
                                    filter_min_total_count = 15, ...) {
  assert_class(design, "FacileDgeModelDefinition")
  if (is.null(prior_count)) prior_count <- 2
  y.all <- as.DGEList(xsamples, assay_name = assay_name, covariates = xsamples)
  # The roughly approximated norm.factors already present in a faciledatastore
  # should be good enough for initial low-pass filtering.
  # y.all <- calcNormFactors(y.all)

  des.matrix <- design(design)[colnames(y.all),,drop=FALSE]
  y.all$design <- des.matrix

  # Remove genes according to `filter` specificaiton
  if (is.character(filter)) {
    if (length(filter) == 1L && filter == "default") {
      # keep only the columns in the design matrix that correspond to the
      # main groups. If this is a ttest, its des.matrix[, design$text_covs].
      # If anova, then we also include the intercept column
      filter.cols <- design$test_covs
      if (is.anova(design)) {
        filter.cols <- c(
          grep("(Intercept)", colnames(des.matrix), ignore.case = TRUE),
          filter.cols)
      }
      keep <- filterByExpr(y.all, des.matrix[, filter.cols],
                           min.count = filter_min_count,
                           min.total.count = filter_min_total_count, ...)
    } else {
      keep <- rownames(y.all) %in% filter
    }
  } else {
    keep <- rep(TRUE, nrow(y.all))
  }

  keep_fraction <- mean(keep)
  keep_n <- sum(keep)

  if (keep_fraction < 0.50 * nrow(y.all)) {
    msg <- glue("Only {format(keep_fraction * 100, digits = 4)}% ",
                "({keep_n}) of features are retained after filtering.")
    warnings <- c(warnings, msg)
  }

  y <- y.all[keep,,keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  if (method == "edgeR-qlf") {
    out <- estimateDisp(y, y[["design"]], robust = TRUE)
  } else if (method == "voom") {
    if (with_sample_weights) {
      out <- voomWithQualityWeights(y, y[["design"]], ...)
    } else {
      out <- voom(y, y$design, save.plot = TRUE, ...)
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

  dropped <- setdiff(rownames(y.all), rownames(out))
  out
}

#' @noRd
#' @importFrom stats hat
#' @importFrom limma arrayWeights
.biocbox_create_EList <- function(xsamples, assay_name, assay_type,
                                  design, filter, method, with_sample_weights,
                                  prior_count = 0.25,
                                  filter_min_expr = 1, ...) {
  if (is.null(prior_count)) prior_count <- 0.25

  # assay_type is one of c("normcounts","lognorm")
  # Note that as.DGEList always assembles the assay matrix with
  # normalized = FALSE
  y.all <- as.DGEList(xsamples, assay_name = assay_name, covariates = xsamples)
  e <- y.all[["counts"]]
  if (assay_type %in% c("normcounts")) {
    e <- log2(e + prior_count)
  }

  # Remove genes according to `filter` specificaiton
  if (is.character(filter)) {
    if (length(filter) == 1L && filter == "default") {
      # if (assay_type == "lognorm") {
      #   min.expr <- log2(2)
      # } else {
      #   min.expr <- 1
      # }
      min.expr <- filter_min_expr
      # Code below taken from edgeR:::filterByExpr.default function
      h <- hat(design)
      min.samples <- 1 / max(h)
      if (min.samples > 10) {
        min.samples <- 10 + (min.samples - 10) * 0.7
      }
      tol <- 1e-14
      keep <- rowSums(e >= min.expr) >= (min.samples - tol)
      # end filterByExpr block
      elist <- elist[keep,]
    } else {
      keep <- rownames(y.all) %in% filter
    }
  } else {
    keep <- rep(TRUE, nrow(y.all))
  }

  elist <- list()
  elist[["E"]] <- e[keep,,drop = FALSE]
  elist[["genes"]] <- y[["genes"]][keep,,drop = FALSE]
  elist[["targets"]] <- y[["samples"]]
  elist <- new("EList", elist)
  elist[["design"]] <- design

  if (with_sample_weights) {
    elist <- arrayWeights(elist, elist[["design"]])
  }

  dropped <- setdiff(rowanmes(y.all), rownames(elilst))

  elist
}
