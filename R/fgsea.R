#' Performs Gene Set Enrichment Analyses
#'
#' For now we only support GSEA methods that can be run from a
#' `FacileDGETestResult`, such as `"cameraPR"` and `"goseq"`
#'
#' @export
#' @importFrom multiGSEA GeneSetDb multiGSEA
#'
#' @param x A `FacileAnalysisResult` object
#' @param gdb A `multiGSEA::GeneSetDb` object
#' @param methods the GSEA methods to use on `x`.
#' @return A FacileGSEAResult object, which includes a MultiGSEAResult object
#'   as it's `result()`.
fgsea <- function(x, gdb, methods, ...) {
  UseMethod("fgsea", x)
}

#' @noRd
#' @export
fgsea.FacileTtestDGEResult <- function(x, gdb, methods = c("cameraPR", "goseq"),
                                       min_logFC = 1, max_padj = 0.10, ...) {

}

#' @noRd
#' @export
fgsea.FacileAnovaDGEResult <- function(x, gdb, methods = "goseq",
                                       max_padj = 0.10, ...) {

}

#' @noRd
#' @export
fgsea.FacilePCAResult <- function(x, methods = "cameraPR", pcs = NULL, ...) {

}

#' @section Accessing Results:
#' A call to `result(fgea.res)` will return the `multiGSEA::MultiGSEAResult`
#' object. If the user specifies the name of a GSEA method used in the analysis,
#' then the summarized results from that method will be returned by passing
#' through to the `multiGSEA::result`function, ie. `r1` and `r2` in the code
#' below are equivalent:
#'
#' ```r
#' r1 <- result(fgsea.res, name = "cameraPR")
#' r2 <- result(fgsea.res) %>% multiGSEA::result("cameraPR")
#' ```
#'
#' @rdname fgsea
#' @export
result.FacileGSEAResult <- function(x, name = "result", ...) {
  mgres <- x[["result"]]
  if (name == "result") {
    out <- mgres
  } else {
    out <- multiGSEA::result(mgres, name, ...)
  }
  out
}
