#' Peforms a differential expression analysis and GSEA.
#'
#' @export
#' @importFrom multiGSEA calculateIndividualLogFC
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

  if (!xor(is.null(coef), is.null(contrast))) {
    msg <- "Only either coef or conrtrast can be non-NULL"
    errors <- c(errors, msg)
  }

  .fds <- assert_class(x$fds, "FacileDataStore")

  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
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
    y <- fdge_biocbox(x, assay_name, method, dge_methods, filter, ...)
    messages <- c(messages, attr(y, "messages"))
    warnings <- c(warnings, attr(y, "warnings"))
    errors <- c(errors, attr(y, "errors"))
  } else {
    y <- NULL
    gsea <- NULL
    dge <- NULL
  }

  method <- y$dge_method
  method_info <- filter(dge_methods, dge_method == method)

  testme <- if (is.null(coef)) contrast else coef
  dge <- calculateIndividualLogFC(y, y$design, contrast = testme)
  gsea <- NULL

  out <- list(
    # Standard FacileAnalysisResult things
    fds = .fds,
    biocbox = y,
    dge = dge,
    gsea = gsea,
    messages = messages,
    warnings = warnings,
    errors = errors)

  clazz <- switch(class(x)[1L],
                  FacileTtestDGEModelDefinition = "FacileTtestResult",
                  FacileAnovaModelDefinition = "FacileAnovaResult",
                  NULL)
  class(out) <- c(clazz, "FacileDGEResult", "FacileAnalysisResult")
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

#' A table of assay_type,dge_method combination parameters
#'
#' This table is used to match the assay_type,dge_method combination with
#' the appropriate bioc container class `bioc_class`, and default edgeR/limma
#' model fitting params.
#'
#' @noRd
fdge_methods <- function(assay_type = NULL) {
  # assay_type values : rnaseq, umi, affymrna, affymirna, log2

  # This is a table of assay_type : dge_method possibilites. The first row
  # for each assay_type is the default analysis method
  assay_methods <- tribble(
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

  method_params <- tribble(
    ~dge_method,  ~robust_fit,  ~robust_ebayes,  ~trend_ebayes,
    "voom",       FALSE,        FALSE,           FALSE,
    "qlf",        TRUE,         FALSE,           FALSE,
    "trended",    FALSE,        FALSE,           TRUE,
    "limma",      FALSE,        FALSE,           FALSE)

  info <- left_join(assay_methods, method_params, by = "dge_method")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}
