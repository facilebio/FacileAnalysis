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
#' @importFrom glue glue
#' @importFrom FacileViz plot
#' @importFrom broom tidy
NULL

# Externally defined generics to re-export on package load ---------------------

#' @importFrom broom tidy
#' @export
broom::tidy

#' @importFrom FacileData %>%
#' @export
FacileData::`%>%`

# dplyr API

# @export distinct
# @export mutate
# @export select
# @export filter
# @export summarise
# @export summarize
# @export arrange
# NULL
