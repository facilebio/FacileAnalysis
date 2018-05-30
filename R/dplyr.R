# immersive results (iresult objects) all have a `$df` element, which
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
#' @method mutate FacileAnalysis
mutate.FacileAnalysis <- function(.data, ...) mutate(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method select FacileAnalysis
select.FacileAnalysis <- function(.data, ...) select(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method filter FacileAnalysis
filter.FacileAnalysis <- function(.data, ...) filter(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method summarise FacileAnalysis
summarise.FacileAnalysis <- function(.data, ...) summarise(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method summarize FacileAnalysis
summarize.FacileAnalysis <- function(.data, ...) summarise(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method arrange FacileAnalysis
arrange.FacileAnalysis <- function(.data, ...) arrange(tidy(.data), ...)

#' @rdname dplyr-api
#' @export
#' @method group_by FacileAnalysis
group_by.FacileAnalysis <- function(.data, ...) group_by(tidy(.data), ...)
