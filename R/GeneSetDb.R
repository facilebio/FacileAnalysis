#' Create GeneSetDb objects from different inputs
#'
#' @export
GeneSetDb <- function(x, ...) {
  UseMethod("GeneSetDb", x)
}

#' @noRd
#' @export
GeneSetDb.default <- function(x, ...) {
  args <- list(...)
  multiGSEA::GeneSetDb(x, featureIdMap = args[["featureIdMap"]],
                       collectionName = args[["collectionName"]], ...)
}

#' @noRd
#' @export
GeneSetDb.FacileFeatureSignature <- function(x, ...) {
  multiGSEA::GeneSetDb(x)
}
