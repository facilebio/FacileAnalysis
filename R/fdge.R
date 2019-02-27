#' Peforms a differential expression analysis and GSEA.
#'
#' @export
#'
#' @param x a data source
#' @param numer character vector defining the covariate/groups that
#'   make up the numerator
#' @param denom character vector defining the covariate/groups that
#'   make up the denominator
#' @param fixed character vector defining the covariate/groups to
#'   use as fixed effects
#' @param covariates a data.frame of covariates (columns) for each
#'   sample (rows) in `x`
fdge <- function(x, design = NULL, contrast = NULL, covariates = NULL,
                filter = NULL, method = .dge_methods,
                gsea = "camera", gdb = NULL, ...) {
  UseMethod("fdge", x)
}

#' @export
fdge.FacileAnovaModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                            gsea = NULL, ...) {
  res <- NextMethod(coef = x$coef)
}

#' @export
fdge.FacileTtestDGEModelDefinition <- function(x, assay_name = NULL,
                                               method = NULL, gsea = "cameraPR",
                                               ...) {
  res <- NextMethod(contrast = x$contrast)
}

fdge.FacileDGEModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                          gsea = NULL, coef = NULL,
                                          contrast = NULL, filter = NULL, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  .fds <- assert_class(x$fds, "FacileDataStore")

  if (is.null(assay_name)) {
    assay_name <- default_assay(x)
  }

  ainfo <- try(assay_info(.fds, assay_name), silent = TRUE)
  if (is(ainfo, "try-error")) {
    msg <- glue("assay_name `{assay_name}` not present in FacileDataStore")
    errors <- c(errors, msg)
    assay_type <- "ERROR"
  } else {
    assay_type <- ainfo$assay_type
  }

  dge_methods <- fdge_methods(assay_type)
  if (nrow(dge_methods) == 0L) {
    msg <- glue("Differential expression method `{method}` for assay_type ",
                "`{assay_type}` not found")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0L) {
    y <- bioc_container(x, assay_name, method, dge_methods, filter, ...)
    messages <- c(messages, attr(y.all, "messages"))
    warnings <- c(warnings, attr(y.all, "warnings"))
    errors <- c(errors, attr(y.all, "errors"))
  } else {
    gsea <- NULL
    dge <- NULL
  }

  out <- list(
    # Standard FacileAnalysisResult things
    fds = .fds,
    y = y,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(out) <- c(clazz, "FacileDGEResult", "FacileAnalysisResult")
  out
}

#' Utility function that returns the bioc expression container for the samples
#' enumerated in `sample_info` with the covariates defined there.
#'
#' TODO: the bioc_container function need to adaptively build the right
#' container given the assay_name and method combinations
#'
#' @noRd
#' @importFrom edgeR filterByExpr calcNormFactors estimateDisp
#' @importFrom limma voom
#' @param sample_info a `facile_frame` that enumerates the samples to fetch
#'   data for, as well as the covariates used in downstream analysis
#' @param assay_name the name of the assay to pull data for
#' @param method the name of the dge method that will be used. This will dictate
#'   the post-processing of the data
bioc_container <- function(mdef, assay_name, method, dge_methods = NULL,
                           filter = NULL, ...) {
  assert_class(mdef, "FacileDGEModelDefinition")
  si <- assert_class(mdef$covariates, "facile_frame")
  .fds <- assert_class(fds(mdef), "FacileDataStore")

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
  y.all <- calcNormFactors(y.all)
  y.all$design <- mdef$design[colnames(y.all),]

  if (is.null(filter)) {
    filter <- filterByExpr(y.all, y.all$design, ...)
  }
  assert_logical(filter, len = nrow(y.all))
  fraction_kept <- mean(filter)
  if (fraction_kept < 0.50 * nrow(y.all)) {
    msg <- glue("Only {format(fraction_kept * 100, digits = 4)}% ",
                "({sum(filter)}) of features are retained ",
                "after filtering.")
    warnings <- c(warnings, msg)
  }

  out <- y.all[filter,,keep.lib.sizes = FALSE]
  out <- calcNormFactors(out)

  if (method %in% c("voom", "trended")) {
    out <- voom(out, out$design)
    if (method == "trended") {
      out$weights <- NULL
    }
  } else if (method == "qlf") {
    out <- estimateDisp(out, out$design, rovbust = TRUE)
  }

  attr(out, "messages") <- messages
  attr(out, "warnings") <- warnings
  attr(out, "errors") <- errors
  out
}

tidy.FacileDGEResult <- function(x, result = "dge") {

}

fdge.matrix <- function(x, design = NULL, coef = NULL, contrast = NULL,
                        covariates = NULL, filter = rep(TRUE, nrow(x)),
                        method = .dge_methods,
                        gsea = "camera", gdb = NULL, ...) {
  stopifnot(
    is.data.frame(covariates),
    nrow(covariates) == ncol(x))
  stopifnot(
    all(numer %in% colnames(covariates)),
    all(denom %in% colnames(covariates)),
  )
}

fdge.DGEGLM <- function(x, design = NULL, coef = NULL, contrast = NULL,
                        covariates = NULL, filter = rep(TRUE, nrow(x)),
                        method = .dge_methods,
                        gsea = "camera", gdb = NULL, ...) {

}

fdge.DGEList <- function(x, design = NULL, coef = NULL, contrast = NULL,
                         covariates = NULL, filter = rep(TRUE, nrow(x)),
                         method = .dge_methods,
                         gsea = "camera", gdb = NULL, ...) {

}

fdge.MArrayLM <- function(x, design = NULL, coef = NULL, contrast = NULL,
                          covariates = NULL, filter = rep(TRUE, nrow(x)),
                          method = .dge_methods,
                          gsea = "camera", gdb = NULL, ...) {

}

fdge.EList <- function(x, design = NULL, coef = NULL, contrast = NULL,
                       covariates = NULL, filter = rep(TRUE, nrow(x)),
                       method = .dge_methods,
                       gsea = "camera", gdb = NULL, ...) {

}

# Helpers ======================================================================

.dge_methods <- c("voom", "qlf", "trended", "limma")

#' A table of assay_type to valid dge methods
#' @noRd
fdge_methods <- function(assay_type = NULL) {
  # assay_type values : rnaseq, umi, affymrna, affymirna, log2

  # This is a table of assay_type : dge_method possibilites. The first row
  # for each assay_type is the default analysis method
  info <- tribble(
    ~assay_type, ~dge_method, ~bioc_class,
    "rnaseq",    "voom",      "DGEList",
    "rnaseq",    "qlf",       "DGEList",
    "rnaseq",    "trended",   "DGEList",
    "umi",       "voom",      "DGEList",
    "umi",       "qlf",       "DGEList",
    "umi",       "trended",   "DGEList",
    "tpm",       "trended",   "EList",
    "affymrna",  "limma",     "EList",
    "affymirna", "limma",     "EList",
    "log2",      "limma",     "EList")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}
