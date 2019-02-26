.dge_methods <- c("voom", "qlf", "trended", "limma")

#' A table of assay_type to valid dge methods
#' @noRd
fdge_methods <- function(assay_type = NULL) {
  # assay_type values : rnaseq, umi, affymrna, affymirna, log2

  # This is a table of assay_type : dge_method possibilites. The first row
  # for each assay_type is the default analysis method
  info <- tribble(
    ~assay_type, ~dge_method,
    "rnaseq",    "voom",
    "rnaseq",    "qlf",
    "rnaseq",    "trended",
    "umi",       "voom",
    "umi",       "qlf",
    "umi",       "trended",
    "tpm",       "trended",
    "affymrna",  "limma",
    "affymirna", "limma",
    "log2",      "limma")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}

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

fdge.FacileDGEModelDefinition <- function(x, ...) {

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
