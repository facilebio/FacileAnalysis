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

#' @section Linear Model Definitions:
#' This function accepts a model defined using using [fdge_model_def()] and
#' creates the appropriate Bioconductor assay container to test the model
#' given the `assay_name` and dge `method` specified by the user.
#'
#' This function only supports RNA-seq like data, ie. data from assays from
#' these `assay_type`s:
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
#' @importFrom edgeR filterByExpr calcNormFactors estimateDisp
#' @importFrom limma arrayWeights voom voomWithQualityWeights
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
#' @param ... passed down to modeling functions where appropriate.
#' @return a DGEList or EList with assay data in the correct place, and all of
#'   the covariates in the `$samples` or `$targerts` data.frame that are requied
#'   to test the model in `mdef`.
biocbox.FacileDGEModelDefinition <- function(x, assay_name = NULL,
                                             method = NULL,
                                             dge_methods = NULL,
                                             filter = "default",
                                             with_sample_weights = FALSE,
                                             prior_count = 3, ...) {
  assert_class(x, "FacileDGEModelDefinition")
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
  if (!ainfo$assay_type %in% c("rnaseq", "umi", "tpm")) {
    errors <- "DGE analysis only implemented for bulk-rnaseq-like data"
    return(out)
  }

  if (test_string(filter) && filter != "default") {
    errors <- c(
      glue("Invalid `filter` value (`{filter}`). The only valid string value ",
            "the `filter` argument is 'default'"))
    return(out)
  }

  if (is.null(dge_methods)) {
    dge_methods <- fdge_methods(ainfo$assay_type)
  }
  if (is.null(method)) {
    method <- dge_methods$dge_method[1L]
  }

  if (!method %in% dge_methods$dge_method) {
    default_method <- dge_methods$dge_method[1L]
    msg <- glue("Requested dge_method `{method}` not found, using ",
                "`{default_method}` instead")
    warnings <- c(warnings, msg)
    method <- default_method
  }

  # Refactor ###################################################################
  # NOTE: The code from here down assumes we are working with a "count"-like
  # assay, and defaulting to using a DGEList. This is where we need to start
  # hacking to support the other types of assay data.

  y.all <- as.DGEList(si, assay_name = assay_name, covariates = si)
  # The roughly approximated norm.factors already present in a faciledatastore
  # should be good enough for initial low-pass filtering.
  # y.all <- calcNormFactors(y.all)
  y.all$design <- x$design[colnames(y.all),]

  # Remove genes according to `filter` specificaiton
  if (is.character(filter)) {
    if (length(filter) == 1L && filter == "default") {
      keep <- filterByExpr(y.all, y.all$design, ...)
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
    out[["biocbox"]] <- estimateDisp(y, y[["design"]], robust = TRUE)
  } else if (method == "voom") {
    if (with_sample_weights) {
      out[["biocbox"]] <- voomWithQualityWeights(y, y[["design"]], ...)
    } else {
      out[["biocbox"]] <- voom(y, y$design, save.plot = TRUE, ...)
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
    out[["biocbox"]] <- elist
  } else {
    stop("How did we get here?")
  }

  out[["dge_method"]] <- method
  out
}

#' @noRd
#' @export
design.BiocBox <- function(x, ...) {
  box <- x[["biocbox"]]
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
