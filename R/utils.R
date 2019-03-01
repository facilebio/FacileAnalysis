# A place for random utility functions

#' @noRd
#' @export
is.categorical <- function(x, ...) {
  is(x, "character") || is(x, "factor")
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
