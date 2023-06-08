# Packages imported completely for internal use --------------------------------

#' @import dplyr
#' @import checkmate
#' @import FacileData
#' @import tidyr
NULL

# Re-export from FacileData ====================================================

#' @export samples
FacileData::samples

#' @export features
FacileData::features

# Ubiquitously used methods from other package =================================

#' @importFrom sparrow failWith
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
