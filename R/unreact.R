#' Replace ReactiveFacileDataStore versions with their "inert" version.
#'
#' When analyses are run through a gadget, a ReactiveFacileDataStore is stored
#' everywhere a FacileDataStore would be (which is to say: in a lot of places).
#' This function tries to dig through the embedded results of an analysis and
#' replaces the ReactiveFaileDataStore objects with their inner FacileDataStore.
#'
#' This function needs to be called within the gadget, ie. I'm pretty sure it
#' needs to fire in a reactive environment, which you get in the
#' `observeEvent(intput$done, { ... })` expression.
#'
#' @export
#' @param x stuff and things
#' @return a version of x with an "inert" FacileDataStore where the
#'   ReactiveFacileDataStore once was
unreact <- function(x, ...) {
  UseMethod("unreact", x)
}

#' @noRd
#' @export
unreact.default <- function(x, ...) {
  x
}

#' @noRd
#' @export
unreact.list <- function(x, ...) {
  out <- lapply(x, unreact)
  class(out) <- class(x)
  out
}

#' @noRd
#' @export
unreact.FacileAnalysisResult <- function(x, ...) {
  unreact.list(x, ...)
}

#' @noRd
#' @export
unreact.facile_frame <- function(x, ...) {
  as_facile_frame(x, fds(fds(x)))
}


#' @noRd
#' @export
unreact.ReactiveFacileDataStore <- function(x, ...) {
  fds(x)
}
