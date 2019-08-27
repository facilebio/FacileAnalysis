# Keep track of metadata regarding different aspects of running ffsea on
# different FacileAnalysisResult types.
#
# These are implemented as S3 generics so that FacileAnalysisXXX packages can
# be defloped that generate FacileAnalysisResultXXX result types, and can define
# how ffsea is run on those, as well as the metadata associated with that.
#
# These metadata are mostly helpful for the shiny widgets to know what to do.
.ffsea_methods <- function(type = NULL) {
  opts <- tibble::tribble(
    ~type,        ~method,
    "preranked",  "cameraPR",
    "preranked",  "fgsea",
    "enrichment", "enrichtest")
  if (!is.null(type)) {
    type. <- assert_choice(type, opts[["type"]])
    opts <- filter(opts, type == type.)
  }
  opts
}

#' Returns the types of featureset enrichment methods available over different
#' types of results.
#'
#' @export
ffsea_methods <- function(x, ...) {
  UseMethod("ffsea_methods", x)
}

#' @noRd
#' @export
ffsea_methods.default <- function(x, ...) {
  stop("method not implemented for class: ", class(x)[1L])
}

#' @noRd
#' @export
ffsea_methods.FacileAnalysisResult <- function(x, ...) {
  filter(.ffsea_methods(), type %in% ffsea_method_types(x))
}

#' Returns a table of type,method combinations that can be used in ffsea
#'
#' This defines the type of gsea methods that we can pass to multiGSEA to
#' perfrom GSEA.
#'
#' @export
#' @return a type,method table

#' Returns the types of GSEA methods available over a FacileAnalysisResult
#'
#' Methods are either based on "preranked" or "enrichment", ie. a rank-based
#' method takes a vector of ordered/ranked features, and an enrichment-based
#' method takes the top N features from a list of features.
#'
#' Examples of "preranked" based GSEA methods include `cameraPR` and `fgsea`,
#' while "enrichment" methods include `goseq` and hypergeometic test.
#'
#' @noRd
#' @param x a FacileAnalysisResult
#' @return a character vector that is a subset of `c("preranked", "enrichment")`
ffsea_method_types <- function(x, ...) {
  UseMethod("ffsea_method_types", x)
}

#' @noRd
#' @export
ffsea_method_types.default <- function(x, ...) {
  stop("ffsea_method_types not defined for: ", class(x)[1L])
}

#' @noRd
#' @export
ffsea_method_types.FacileTtestAnalysisResult <- function(x, ...) {
  c("preranked", "enrichment")
}

#' @noRd
#' @export
ffsea_method_types.FacileAnovaAnalysisResult <- function(x, ...) {
  c("enrichment")
}

#' @noRd
#' @export
ffsea_method_types.FacilePcaAnalysisResult <- function(x, ...) {
  c("preranked")
}

#' Returns default arguments used to run ffsea over a FacileAnalysisResult
#'
#' Note: these implementation can probably be implemented with a smart use of
#' `formals(x)`, but ...
#'
#' @noRd
#' @export
#' @param x A FacileAnalysisResult
default_ffsea_args <- function(x, ...) {
  UseMethod("default_ffsea_args", x)
}

#' @noRd
#' @export
default_ffsea_args.default <- function(x, ...) {
  stop("default_ffea_args method not defined for: ", class(x)[1L])
}

#' @noRd
#' @export
default_ffsea_args.FacileTtestAnalysisResult <- function(x, ...) {
  list(
    min_logFC = 1, max_padj = 0.10,
    rank_by = "logFC", signed = TRUE)
}

#' @noRd
#' @export
default_ffsea_args.FacileAnovaAnalysisResult <- function(x, ...) {
  list(max_padj = 0.10)
}

#' @noRd
#' @export
default_ffsea_args.FacilePcaAnalysisResult <- function(x, ...) {
  list(dim = 1, signed = TRUE)
}
