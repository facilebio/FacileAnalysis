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
#' @importFrom limma voom
#' @param sample_info a `facile_frame` that enumerates the samples to fetch
#'   data for, as well as the covariates used in downstream analysis
#' @param assay_name the name of the assay to pull data for
#' @param method the name of the dge method that will be used. This will dictate
#'   the post-processing of the data
#' @param filter A filtering policy to remove unintereesting genes.
#'   If `"default"` (which is the default), then [edgeR::filterByExpr()] is
#'   used if we are materializing a `DGEList`, otherwise lowly expressed
#'   features are removed. If `FALSE`, then no filtering is done.
#' @return a DGEList or EList with assay data in the correct place, and all of
#'   the covariates in the `$samples` or `$targerts` data.frame that are requied
#'   to test the model in `mdef`.
biocbox.FacileDGEModelDefinition <- function(x, assay_name, method,
                                             dge_methods = NULL,
                                             filter = "default",
                                             prior_count = 1, ...) {
  assert_class(x, "FacileDGEModelDefinition")
  si <- assert_class(x$covariates, "facile_frame")
  .fds <- assert_class(fds(x), "FacileDataStore")

  messages <- character()
  warnings <- character()
  errors <- character()

  ainfo <- assay_info(.fds, assay_name)
  if (!ainfo$assay_type %in% c("rnaseq", "umi", "tpm")) {
    # TODO: We can implement this for microarrays very easily
    stop("DGE analysis not yet implemented for bulk-rnaseq-like data")
  }

  if (is.null(dge_methods)) {
    dge_methods <- fdge_methods(ainfo$assay_type)
  }

  if (!method %in% dge_methods$dge_method) {
    default_method <- dge_methods$dge_method[1L]
    msg <- glue("Requested dge_method `{method}` not found, using ",
                "`{default_method}` instead")
    warnings <- c(warnings, msg)
    method <- default_method
  }

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

  out <- y.all[keep,,keep.lib.sizes = FALSE]
  out <- calcNormFactors(out)

  if (method %in% c("voom", "trended")) {
    out <- voom(out, out$design)
    if (method == "trended") {
      out$weights <- NULL
    }
  } else if (method == "qlf") {
    out <- estimateDisp(out, out$design, robust = TRUE)
  }

  attr(out, "messages") <- messages
  attr(out, "warnings") <- warnings
  attr(out, "errors") <- errors
  out$dge_method <- method
  out
}
