# A place for random utility functions

#' Extract a feature_id character vector from a feature descriptor, or NULL
#'
#' All of the functions that take in a `features` parameter can call this
#' function to do the id extraction ... I wrote this code a lot and all
#' over the place.
#'
#' @noRd
#' @param x NULL, character string, or tibble with feature_id column
#' @return a character vector of feature ids
extract_feature_id <- function(x, as_tibble = FALSE, ...) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- x[["feature_id"]]
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) {
    stop("Could not extract feature_id from feature descriptor", call. = FALSE)
  }
  if (as_tibble) {
    x <- tibble(feature_id = x)
  }
  x
}

#' Convenience wrapper to require specified packages
#'
#' @noRd
#' @param pkg A character vector of packages to require
#' @param quietly defaults to true
#' @param ... passed into [requireNamespace()]
reqpkg <- function(pkg, quietly = TRUE, ...) {
  assert_character(pkg)
  for (p in pkg) {
    if (!requireNamespace(p, ..., quietly = quietly)) {
      stop("'", p, "' package required, please install it.", call. = FALSE)
    }
  }
}
