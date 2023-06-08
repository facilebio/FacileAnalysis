#' @noRd
is.ttest <- function(x, ...) {
  if (is(x, "ReactiveFacileAnalysisResult")) {
    x <- faro(x)
  }
  is(x, "FacileTtestAnalysisResult") || is(x ,"FacileTtestModelDefinition")
}

#' @noRd
is.anova <- function(x, ...) {
  if (is(x, "ReactiveFacileAnalysisResult")) {
    x <- faro(x)
  }
  is(x, "FacileAnovaAnalysisResult") || is(x ,"FacileAnovaModelDefinition")
}
