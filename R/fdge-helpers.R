#' Checks to see if fdge analysis results are from a ttest or anova
#' 
#' @rdname fdge-helpers
#' @export
is_ttest <- function(x, ...) {
  UseMethod("is_ttest", x)
}

#' @noRd
#' @export
is_ttest.default <- function(x, ...) {
  is(x, "FacileTtestAnalysisResult") || is(x ,"FacileTtestModelDefinition")
}


#' @rdname fdge-helpers
#' @export
is_anova <- function(x, ...) {
  UseMethod("is_anova", x)
}

#' @noRd
#' @export
is_anova.default <- function(x, ...) {
  is(x, "FacileAnovaAnalysisResult") || is(x ,"FacileAnovaModelDefinition")
}
