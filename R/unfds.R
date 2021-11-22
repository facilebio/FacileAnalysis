#' Strip and restore FacileDataStore pointers from/to FacileAnalysis Results
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' These functions are here to support serialization. `unfds` removes the link
#' to a FacileDataStore from any object, and `refds` restores it to only the
#' "blessed" objects.
#'
#' These functions shouldn't really be used by users, but I'm having trouble
#' testing a non-exported S3 function in unit tests ...
#'
#' @export
#' @param x An object to remove the FacileDataStore attribute from
#' @return a stripped out object
unfds <- function(x, ...) {
  UseMethod("unfds", x)
}

#' @export
#' @noRd
unfds.default <- function(x, ...) {
  lifecycle::signal_stage("experimental", "unfds()")
  set_fds(x, NULL)
}

#' @export
#' @noRd
unfds.list <- function(x, ...) {
  classes <- class(x)
  x <- lapply(x, unfds, ...)
  class(x) <- unique(classes)
  set_fds(x, NULL)
}

#' @export
#' @noRd
#' @method unfds data.frame
unfds.data.frame <- function(x, ...) {
  set_fds(x, NULL)
}

#' @export
#' @noRd
unfds.tbl <- function(x, ...) {
  set_fds(x, NULL)
}

#' @export
#' @noRd
unfds.FacileAnalysisResult <- function(x, ...) {
  x[["fds"]] <- NULL
  out <- unfds.list(x, ...)
  out
}

# refds ------------------------------------------------------------------------

#' Reconstitute FacileDataStore links from unfds'd objects
#'
#' Not exported on purpose
#'
#' @export
#' @param an object to reconstitute fds links to
#' @param fds a FacileDataStore
refds <- function(x, fds, ...) {
  UseMethod("refds", x)
}

#' @export
#' @noRd
refds.default <- function(x, fds, ...) {
  x
}

#' @export
#' @noRd
refds.list <- function(x, fds, ...) {
  classes <- class(x)
  out <- lapply(x, refds, fds, ...)
  class(out) <- unique(classes)
  out
}

#' @export
#' @noRd
refds.FacileAnalysisResult <- function(x, fds, ...) {
  assert_class(fds, "FacileDataStore")
  x <- refds.list(x, fds, ...)
  x[["fds"]] <- fds
  x
}

#' @export
#' @noRd
refds.facile_frame <- function(x, fds, ...) {
  set_fds(x, fds)
}
