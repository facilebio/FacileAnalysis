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

#' @noRd
#' @export
#' @importFrom FacileShine initialized
initialized.ReactiveGeneSetDb <- function(x, ...) {
  is(x$gdb(), "GeneSetDb") && nrow(x$geneSets()) > 0
}

#' Returns the filtered unreactive GeneSetDb from the reactiveGeneSetDb module.
#'
#' The user can choose to include subsets of collections, as well as
#' genesets of a certain size.
#'
#' This can only be run within a reactive context.
#'
#' @noRd
#' @export
#' @importFrom multiGSEA geneSets
GeneSetDb.ReactiveGeneSetDb <- function(x, ...) {
  gdb <- x$gdb()
  gsets.all <- geneSets(gdb)
  gsets.selected <- x$geneSets()

  keep <- is.element(
    paste(gsets.all$collection, gsets.all$name),
    paste(gsets.selected$collection, gsets.selected$name))

  if (!all(keep)) {
    gdb <- gdb[keep]
  }

  gdb
}
