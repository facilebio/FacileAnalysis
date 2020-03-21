# A place for random utility functions

#' Extract a feature_id characger vector from a feature descriptor, or NULL
#'
#' All of the functions that take in a `features` parameter can call this
#' function to do the id extraction ... I wrote this code a lot and all
#' over the place.
#'
#' @noRd
#' @param x NULL, character string, or tibble with feature_id column
extract_feature_id <- function(x, ...) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- x[["feature_id"]]
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) {
    stop("Could not extract feature_id from feature descriptor", call. = FALSE)
  }
  x
}

#' A wrapper to bscols that suppresswarnings.
#'
#' I am using bscols to create multiple rows, and in this case it always
#' throws a warning, which I want to suppress
#'
#' @noRd
#' @importFrom crosstalk bscols
bscols. <- function(..., suppress_warning = TRUE) {
  wrap. <- if (suppress_warning) suppressWarnings else identity
  wrap.(bscols(...))
}
