# FacileAnalysisResult objects =================================================

#' @export
#' @rdname getsetdb
fds.FacileAnalysisResult <- function(x) {
  return(x[["fds"]])
}

# dplyr (FacileAnalysisResult) -------------------------------------------------
#
# FacileAnalysisReult have data to extract, they should be in tidy forma
# is a data.frame, therefore the results should be chainable to a dplyr function

# According to the documentation (http://dplyr.tidyverse.org/) the dplyr
# API is defined by these functions, so let's see if we can do something
# useful with them
#
# * mutate(): adds new variables that are functions of existing variables
# * select(): picks variables based on their names.
# * filter(): picks cases based on their values.
# * summarise(): reduces multiple values down to a single summary.
# * arrange(): changes the ordering of the rows.

#' ImmersiveAnalysis result objects can be maniuplated with dplyr verbs
#'
#' @rdname dplyr-api
#' @export
#' @method mutate FacileAnalysisResult
mutate.FacileAnalysisResult <- function(.data, ...) mutate(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method select FacileAnalysisResult
select.FacileAnalysisResult <- function(.data, ...) select(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method filter FacileAnalysisResult
filter.FacileAnalysisResult <- function(.data, ...) filter(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method summarise FacileAnalysisResult
summarise.FacileAnalysisResult <- function(.data, ...) summarise(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method summarize FacileAnalysisResult
summarize.FacileAnalysisResult <- function(.data, ...) summarise(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method arrange FacileAnalysisResult
arrange.FacileAnalysisResult <- function(.data, ...) arrange(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method group_by FacileAnalysisResult
group_by.FacileAnalysisResult <- function(.data, ...) group_by(tidy(.data), ...)

# broom (FacileAnalysisResult) -------------------------------------------------

#' Extract a tidy data.frame from an `iresult`
#'
#' @export
#' @method tidy FacileAnalysisResult
#'
#' @param x an immersive result object (`iresult`)
#' @return a tidy data.frame of the results from an immersive analysis
tidy.FacileAnalysisResult <- function(x, ...) x$tidy

