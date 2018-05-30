#' Launch a shiny gadget to explore an immersive result
fxplore <- function(x, ...) UseMethod("fxplore")

#' Creates an interactive plot of the immersive result
fplot <- function(x, ...) UseMethod("fplot")

#' Extract a tidy data.frame from an `iresult`
#'
#' @export
#' @method tidy FacileAnalysis
#'
#' @param x an immersive result object (`iresult`)
#' @return a tidy data.frame of the results from an immersive analysis
tidy.FacileAnalysis <- function(x, ...) x$tidy

