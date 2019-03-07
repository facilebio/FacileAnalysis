#' Takes care of the tedium of wrapping a facilemodule to run as a gadget.
#'
#' @param x A FacileDataStore
#' @param samples an optional subset of samples to launch the gadget over
#' @param module the to launch as a gadget
#' @param ui An optional UI function for the gadget. If none is provided,
#'  it will default to calling `moduleUI()`
#'
#'  @examples
#'  \dontrun{
#'  fds <- FacileData::exampleFacileDataSet()
#'  stuff <- fgadget(fds, fdgeAnalysis)
#'  }
fgadget <- function(x, samples = NULL, module, ui = NULL, ...) {

}
