#' Gadget to run both DGE and GSEA
#'
#' @export
fdgeseaGadget <- function(x, title = "DGE and GSEA",
                          height = 800, width = 1000, ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fDgeSeaAnalysis, fDgeSeaAnalysisUI, x, title = title,
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
#' This also gives us the opportunity to exercise different methods of
#' interaction between independant analysis results. In the DGE and GSEA
#' scenario, for instance, we might want the linked brushing that happens within
#' the `fdge` volcano to do a "significant gene set search" in the fgsea result.
#'
#' TODO: Enable / Disable FSEA widgets and results when DGE result isn't ready
#' or a new one is generated.
#'
#' @export
fDgeSeaAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  model <- callModule(fdgeModelDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)

  # Something to think about given the current design of these modules is if we
  # can broadcast the feature brushing that happens here in such a way that
  # a "listner" in the ffseaView model can react to it. Likewise can the
  # ffseaView module broadcast feature brushing in such a way that the volcano
  # module in the fdgeView module can react to.
  dge_view <- callModule(fdgeView, "dge_view", rfds, dge,  ...,
                        feature_selection = session$ns("volcano"),
                        sample_selection = session$ns("samples"),
                        debug = debug)

  fsea <- callModule(ffseaRun, "fsea", rfds, dge, gdb, ..., debug = debug)
  fsea_view <- callModule(ffseaView, "fsea_view", rfds, fsea, ...,
                          debug = FALSE)

  # Only show the view UI when there is (at least) a DGE result.
  observe({
    res. <- req(faro(dge))
    show <- is(res., "FacileDgeAnalysisResult")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    main = list(dge = dge, fsea = fsea),
    .ns = session$ns)

  class(vals) <- c(
    "ReactiveFacileMultiAnalysisResultContainer" # are you being serious?
  )

  vals
}


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
           fdgeModelDefRunUI(ns("model"), debug = debug))),

    # DGE and GSEA parameters, side-by-side
    tags$div(
      id = ns("configbox"),
      fluidRow(
        box.(title = "Differential Expression Parameters", width = 6,
             fdgeRunUI(ns("dge"), debug = debug)),
        box.(title = "Feature Set Enrichment Parameters", width = 6,
             ffseaRunUI(ns("fsea"), debug = debug)))),
    # Result Container
    hidden(
      tags$div(
        id = ns("viewbox"),
        tabBox(
          title = "Analysis Results", width = 12,
          tabPanel("Differential Expression", fdgeViewUI(ns("dge_view"))),
          tabPanel("Feature Set Enrichment", ffseaViewUI(ns("fsea_view")))))
    )
  )
}
