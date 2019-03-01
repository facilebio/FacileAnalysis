#' Performs Gene Set Enrichment Analyses
#'
#' For now we only support GSEA methods that can be run from a
#' `FacileDGETestResult`, such as `"cameraPR"` and `"goseq"`
#'
#' @export
#' @importFrom multiGSEA GeneSetDb multiGSEA
fgsea <- function(x, ..) {
  UseMethod("fgsea", x)
}

fgsea.FacileTtestDGEResult <- function(x, methods = c("cameraPR", "goseq"),
                                       min_logFC = 1, max_padj = 0.10, ...) {

}

fgsea.FacileAnovaDGEResult <- function(x, methods = "goseq", max_padj = 0.10,
                                       ...) {

}
