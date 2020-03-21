#' @include fdge.R
NULL

#' Gadget to run both DGE and GSEA
#'
#' @export
#' @examples
#' \dontrun{
#' gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
#' efds <- FacileData::exampleFacileDataSet()
#' xs <- FacileData::filter_samples(efds, indication == "CRC")
#' fdgeseaGadget(xs, gdb)
#' }
fdgeseaGadget <- function(x, gdb = NULL, title = "DGE and GSEA",
                          height = 800, width = 1000, ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fDgeSeaAnalysis, fDgeSeaAnalysisUI, x, gdb = gdb, title = title,
             height = height, width = width, ...)
}

#' An analysis module that combines differential expression with GESA.
#'
#' This module enables the interactive configuration and execution of
#' "paired" differential expression and feature set enrichment.
#'
#' Although the focus of the analysis workflow in this package is based on
#' singular units of analysis, there are times like this where we almost always
#' want to perform two anlayses together, differential expression and GSEA
#' on the same contrast. This shiny-module presents an interface to both.
#'
#' @section Development Thougts:
#' This also gives us the opportunity to exercise different methods of
#' interaction between independant analysis results. In the DGE and GSEA
#' scenario, for instance, we might want the linked brushing that happens within
#' the `fdge` volcano to do a "significant gene set search" in the fgsea result.
#'
#' More specifically, maybe the volcano plot in the `dge_view` module should
#' be able to broadcast the feature brushing in such a way that the "listener"
#' in the `ffsea_vew` module can react to it. Similarly for the `ffsea_vew`
#' module should be able to broadcast the features that are brushed within
#' it to the `dge_vew` volcano and statistics tables ... or not?
#'
#' @export
fDgeSeaAnalysis <- function(input, output, session, rfds, gdb = NULL, ...,
                            debug = FALSE) {
  # fdge bits ..................................................................
  model <- callModule(flmDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)
  dge_view <- callModule(fdgeView, "dge_view", rfds, dge,  ...,
                        feature_selection = session$ns("volcano"),
                        sample_selection = session$ns("samples"),
                        debug = debug)

  # ffsea bits .................................................................
  fsea <- callModule(ffseaRun, "fsea", rfds, dge, gdb = gdb, ..., debug = debug)
  fsea_view <- callModule(ffseaView, "fsea_view", rfds, fsea, ...,
                          debug = FALSE)

  # toggle UI logic ............................................................
  # Only show the view UI when there is (at least) a DGE result.
  observe({
    res. <- req(faro(dge))
    show <- is(res., "FacileDgeAnalysisResult")
    toggleElement("viewbox", condition = show)
  })

  # TODO: Enable / Disable FSEA widgets and results when DGE result isn't ready
  # or a new one is generated.

  # wrap up ....................................................................
  vals <- list(
    main = list(dge = dge, fsea = fsea),
    .ns = session$ns)

  class(vals) <- c(
    # This is a serious whopper of a name -- are you being serious?
    "ReactiveFacileMultiAnalysisResultContainer"
  )

  vals
}


#' @export
#' @noRd
#' @importFrom shinydashboard tabBox
#' @importFrom shinyjs hidden
#' @importFrom shiny tabPanel
fDgeSeaAnalysisUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  box. <- shinydashboard::box
  tagList(
    # Linear Model Definitions
    tags$div(
      id = ns("modelbox"),
      box.(title = "Model Definition", width = 12,
           flmDefRunUI(ns("model"), debug = debug))),

    # DGE and GSEA parameters, side-by-side
    tags$div(
      id = ns("configbox"),
      fluidRow(
        box.(title = "Differential Expression Parameters", width = 6,
             fdgeRunUI(ns("dge"), debug = debug)),
        box.(title = "Set Enrichment Parameters", width = 6,
             ffseaRunUI(ns("fsea"), debug = debug)))),
    # Result Container
    hidden(
      tags$div(
        id = ns("viewbox"),
        tabBox(
          title = "Analysis Results", width = 12,
          tabPanel("Differential Expression", fdgeViewUI(ns("dge_view"))),
          tabPanel("Set Enrichment", ffseaViewUI(ns("fsea_view")))))
    )
  )
}
