# FacileAnalysis S3 Generics ===================================================

#' Creates an interactive plot of the FacileAnalysisResult
render <- function(x, ...) UseMethod("render", x)

#' Launch a shiny gadget to explore a FacileAnalysisResult
shine <- function(x, ...) UseMethod("shine", x)

# Getters and Setters ----------------------------------------------------------
#' Retrieve the input object from a FacilePlot
#' @export
input_data <- function(x, ...) UseMethod("input_data", x)

#' Extracts the data used for plotting
#' @export
plot_data <- function(x, ...) UseMethod("plot_data", x)

#' Extract a FacileAnalysisResult from something
#' @export
analysis_result <- function(x, ...) UseMethod("analysis_result", x)


# FacileAnalysisResult objects =================================================

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


# Plots ========================================================================

# A FacilePlot object (fP) includes the following elements:
#
# * `fp$facile_analysis`: The `FacileAnalysisResult` object that the plot was
#   generated from. This is made accessible by alling `fanalysis_result(fp)`.
# * `fp$input_data`: A slimmed down version of the data used for the final plot.
#   A bigger version of this data should be available via:
#   `fp %>% fanalysis_result() %>% tidy()`
# * `fp$plot`: The plotly (for now) object that is the actual plot itself. This
#   is made accessible via `plot(fp)` method
#   (TODO: should the plot,FacilePlot method really be called fplot?)
# * `fp$params`: The parameters used to configure the plot, made accessible via
#   `fparams(fp)` (TODO: make fparams S3)

#' Extract (plotly)-plot object from FacilPlots objects
#'
#' This creates an
#'
#' @noRd
#'
#' @export
#' @method plot FacilePlot
plot.FacilePlot <- function(x, y, ...) {
  assert_class(x, "FacilePlot")
  x$plot
}

#' Print a FacilePlot to console/viewer
#' @noRd
print.FacilePlot <- function(x, ..., view = interactive()) {
  # delegates to htmlwidgets::print.htmlwdiget
  invisible(print(plot(x), ..., view = view))
}

#' @export
#' @method analysis_result FacilePlot
analysis_result.FacilePlot <- function(x, ...) {
  assert_class(x, "FacilePlot")
  x$facile_analysis
}

# Letting the knit_print method pass through from
# plotly -> htmlwidgets::knit_print seems to be doing the right thing
# for now
# knit_print.FacilePlot <- function(x, ..., options = NULL) {
#   # delegates to htmlwidgets::knit_print.htmlwdiget
#   knitr::knit_print(plot(x))
# }
#

