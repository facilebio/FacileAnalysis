#' Extracts messages, warnings, and error messages from a FacileAnalysis Result
#' 
#' @export
#' @rdname messages-warnings-errors
messages <- function(x, ...) {
  if (missing(x) || !is(x, "FacileAnalysisResult")) {
    stop("Undefined call to errors")
  }
  x[["messages"]]
}

#' @export
#' @rdname messages-warnings-errors
warnings <- function(...) {
  if (...length() == 0L || !is(..1, "FacileAnalysisResult")) {
    base::warnings(...)
  } else {
    ..1[["warnings"]]
  }
}

#' @export
#' @rdname messages-warnings-errors
errors <- function(x, ...) {
  if (missing(x) || !is(x, "FacileAnalysisResult")) {
    stop("Undefined call to errors")
  }
  x[["errors"]]
}


