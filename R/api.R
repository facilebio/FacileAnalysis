# Universal FacileAnalysisResult methods =======================================

#' Compares two (or more(?)) FacileAnalysisResults against each other.
#'
#' This was initially motivated by the desire to compare the results of two
#' differential expressio analyses against each other. Extending the idea
#' to comparing GSEA results (using Thomas' enrichmentmap idea) is a natural
#' extension. There may be othres.
#'
#' Not every FacileAnalysisResult can be compared, I guess, in which case it
#' we will throw an error
compare <- function(x, y, ...) {
  UseMethod("compare", x)
}

compare.FacileAnalysisResult <- function(x, y, ...) {
  msg <- glue("No `compare` method defined for `{class(x)[1L]}` result type.")
}

# Retrieves the main (or alternative) result from a FacileAnalysis
#
# This should always return a tidy data frame of one sort or another.
#
# FacileAnalysisResults should should all aslo have a top-level result
# accessible via obj[["result"]]
#
#' @importFrom multiGSEA result
#' @export result
NULL

# This is defined in multiGSEA@develop
# result <- function(x, name = "result", ...) {
#   UseMethod("result", x)
# }

#' @export
#' @noRd
result.FacileAnalysisResult <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' Extract the value of a parameter used in a FacileAnalysis result
#'
#' @export
#' @param x A FacileAnalysisResult
#' @param name the name of the parameter to extract
param <- function(x, name, ...) {
  UseMethod("param", x)
}

#' @noRd
#' @export
param.FacileAnalysisResult <- function(x, name, ...) {
  params. <- assert_list(x[["params"]], names = "unique")
  assert_choice(name, names(params.))
  params.[[name]]
}

#' Extracts feature or sample ranks from a FacileAnalysis
#'
#' Among other thigns, analyses can often provide rankings over features or
#' samples. For instance, a differential expression result may provide ranking
#' over the genes. This generic extract ranks of one type or another.
#'
#' @export
#' @param x A `FacileAnalysisResult`
ranks <- function(x, ...) {
  UseMethod("ranks", x)
}

#' Converts an analysis result into a "signature" type of thing.
#'
#' This is mainly for feature rankings, but could be something else? Returns
#' a table of features that can be piped into GeneSetDb() constructor, should
#' insipiration strike.
#'
#' Most of the `signature.*` functions defined over the result of analysis first
#' passes through a call to `ranks()` first, which then calls the `signature()`
#' method on those ranks, ie. the two segments below will produce the same
#' result.
#'
#' ```r
#' sigs1 <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   signature()
#' sigs2 <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   ranks() %>%
#'   signature()
#' ```
#'
#' @export
#' @param x A `FacileAnalysisResult`.
#' @examples
#' pca.sigs <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   signature()
signature <- function(x, ...) {
  UseMethod("signature", x)
}

# Facile API Implementations ===================================================

#' FacileAnalysisResult objects should be able to fetch the FacileDataStore
#' they were created from.
#'
#' @export
#' @noRd
fds.FacileAnalysisResult <- function(x) {
  return(x[["fds"]])
}

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
viz <- function(x, ...) {
  UseMethod("viz", x)
}

#' @noRd
#' @export
viz.FacileGadgetResult <- function(x, ...) {
  viz(result(x), ...)
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

#' @noRd
#' @export
report.FacileGadgetResult <- function(x, ...) {
  report(result(x), ...)
}

# Statistical Modeling Stuff ===================================================

#' Extract the design matrix from objects used in a linear modeling step
#'
#' @export
design <- function(x, ...) {
  UseMethod("design", x)
}

#' Extract the model used for an analysis
#'
#' For a DGE analysis, this is the fdge_model_def object
#'
#' @export
model <- function(x, ...) {
  UseMethod("model", x)
}
