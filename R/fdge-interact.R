# Interactivity and Vizualization over FacileDGEResults ========================

#' @noRd
#' @export
vizualize.FacileTtestDGEResult <- function(x, result = c("features", "gsea"),
                                           with_boxplot = TRUE, ntop = 100,
                                           max_padj = 0.01, min_abs_logFC = 1,
                                           feature = NULL, event_source = "A",
                                           ...) {
  if (FALSE) {
    result <- "features"
    with_boxplot <- TRUE
    topn <- 100
    max_padj <- 0.01
    min_abs_logFC <- 1
  }
  result <- match.arg(result)
  type <- match.arg(type)

}

#' @section Interacting with results:
#'
#' The `report` function will create an htmlwidget which can be explored by
#' the analyst or dropped into an Rmarkdown report.
#'
#' `report(result, "dge", max_padj = 0.05, min_abs_logFC = 1)` will create a
#' side-by-side volcano and datatable for differential expression results.
#'
#' @export
#' @rdname fdge
report.FacileTtestDGEResult <- function(x, result = "dge",
                                        max_padj = 0.01, min_abs_logFC = 1,
                                        event_source = "A", ...) {
  # result <- match.arg(result, c("dge", multiGSEA::resultNames(result$gsea)))
  result <- match.arg(result, "dge")

}
