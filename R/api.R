# Vizualization and Rmarkdon reporting =========================================

#' Methods to interactively explore and report FacileAnalysisResults.
#'
#' The `vizualize`, `shine`, and `report` triumverate provide the analyst with
#' the tools required to interact and explore the results of a FacileAnalysis.
#'
#' @section Vizualize:
#' The `vizualize` functions generate an analysis-specific interactive
#' htmlwidget for the analysist to explore.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#'
#' @param x A `FacileAnalysisResult` object
#' @param ... passed down to the `x`-specific `vizualize.*`, `report.*`, and
#'   `shine.*` functions.
vizualize <- function(x, ...) {
  UseMethod("vizualize", x)
}

#' @section Shine:
#' The `shine` functions generate shiny gadgets that provide a more interactive
#' view over a `FacileAnalysisResult`. This empowers the analyst to provide more
#' context around the results, likely by leveraging all of the data available
#' within the FacileDataStore.
#'
#' The respective `shine` functions must return a `FacileAnalysisShine` object
#' invisibly to the caller. These should be able to be past into an overladed
#' `report` function, the result of which can be embedded into an Rmarkdown
#' report. In this way the analyst can embed a feature-reduced version of what
#' was observed in the gadget into an Rmarkdown report.
#'
#' I'm not sure how exactly we can do this, but perhaps this will require some
#' code generation that the analyst can copy and paste paste into the Rmarkdown
#' document. This might simply be a parameterized version of the `report`
#' function call.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#' @aliases shine
shine <- function(x, ...) {
  UseMethod("shine", x)
}

#' @section Report:
#' The `report` function produces an object that can be embedded into an
#' Rmarkdown document. The implementation of these functions will likely
#' result in a parameterized call to the respective `vizualize` function.
#'
#' Two `report` functions should be created per FacileAnalysis. One that accepts
#' the `FacileAnalysisResult` object itself, and another that accepts the
#' analysis-specific `FacileAnalysisShine` object, which should parameterize
#' the final `visualize` call so that it can be embedded seamlessly into
#' an Rmarkdown report for "offline" viewing.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#' @aliases report
report <- function(x, ...) {
  UseMethod("report", x)
}


# FacileAnalysisResult objects =================================================

#' @export
#' @noRd
fds.FacileAnalysisResult <- function(x) {
  return(x[["fds"]])
}

# broom (FacileAnalysisResult) -------------------------------------------------

#' Extract a tidy data.frame from an `iresult`
#'
#' @export
#' @method tidy FacileAnalysisResult
#'
#' @param x an immersive result object (`iresult`)
#' @return a tidy data.frame of the results from an immersive analysis
tidy.FacileAnalysisResult <- function(x, ...) x[["tidy"]]

