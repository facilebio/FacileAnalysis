#' Adds gadget-specific class info to the class() of a FacileAnalysisResut
#'
#' @noRd
#' @return a character vector of classes for `x`
classify_as_gadget <- function(x, ...) {
  if (is.list(x) && length(x) == 0L) {
    out <- "FailedGadgetResult" # ugh, not really but ...
  } else if (is(x, "FacileAnalysisResult")) {
    oc <- class(x)
    out <- c(
      sub("Result$", "GadgetResult", oc[1L]),
      head(oc, -1L),
      "FacileAnalysisGadgetResult",
      tail(oc, 1L))
    out <- unique(out)
  } else {
    out <- class(x)
  }
  out
}
