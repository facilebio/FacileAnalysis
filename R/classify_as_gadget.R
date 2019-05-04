#' Adds gadget-specific class info to the class() of a FacileAnalysisResut
#'
#' @noRd
classify_as_gadget <- function(x, ...) {
  if (is.list(x) && length(x) == 0L) return("FailedGadgetResult")
  stopifnot(is(x, "FacileAnalysisResult"))
  oc <- class(x)
  out <- c(
    sub("Result$", "GadgetResult", oc[1L]),
    head(oc, -1L),
    "FacileAnalysisGadgetResult",
    tail(oc, 1L))
  out
}
