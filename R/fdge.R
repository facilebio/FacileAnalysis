#' Peforms a differential expression analysis and GSEA.
#'
#'
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

fdge.FacileDataStore <- function(x, design = NULL, coef = NULL, contrast = NULL,
                                 covariates = NULL, filter = rep(TRUE, nrow(x)),
                                 method = .dge_methods,
                                 gsea = "camera", gdb = NULL,
                                 assay_name = default_assay(x), ...) {
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
.dge_methods <- c("qlf", "voom", "trended", "limma")
