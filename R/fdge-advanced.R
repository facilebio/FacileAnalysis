# Manual analyst intervention to run any type of model

#' FacileAnalysis wrapper for a manually performed differential expression.'
#'
#' This allows an analyst to perform a differential expression analysis
#' outside of the facileverse, and wrap it into an object that can be used
#' with the viz,report,shine functions.
as_fdge <- function(x, ...) {

}

if (FALSE) {
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
}
