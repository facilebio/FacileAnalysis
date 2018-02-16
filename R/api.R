#' Launch a shiny gadget to explore an immersive result
ixplore <- function(x, ...) UseMethod("ixplore")

#' Creates an interactive plot of the immersive result
iplot <- function(x, ...) UseMethod("iplot")

#' Extract a tidy data.frame from an `iresult`
#'
#' @export
#' @method tidy iresult
#'
#' @param x an immersive result object (`iresult`)
#' @return a tidy data.frame of the results from an immersive analysis
tidy.iresult <- function(x, ...) x$tidy

