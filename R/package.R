# Packages imported completely for internal use --------------------------------

#' @import dplyr
#' @import checkmate
#' @import FacileData
#' @import tidyr
NULL

#' @export samples

# Ubiquitously used methods from other package =================================

#' @importFrom FacileShine annotation
#' @export annotation
NULL

#' @importFrom multiGSEA result
#' @export result
multiGSEA::result

#' @importFrom glue glue
#' @importFrom FacileViz plot
NULL

#' @importFrom FacileViz plot_data
#' @export
FacileViz::plot_data

# Externally defined generics to re-export on package load ---------------------

#' @importFrom broom tidy
#' @export
broom::tidy

#' @importFrom FacileData %>%
#' @export
FacileData::`%>%`

#' @importFrom shiny callModule
#' @export
shiny::callModule
