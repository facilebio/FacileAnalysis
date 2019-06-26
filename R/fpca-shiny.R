#' Perform a PCA analysis on a (subset) of a FacileDataStore
#'
#' Interactively runs PCA on a FacileDataStore. If you want to run PCA on a
#' subset of samples, pass in a facile_frame sample descriptor
#'
#' @export
#' @examples
#' if (interactive()) {
#' efds <- exampleFacileDataSet()
#' # run tumor vs normal comparisons vs each, then run compare9) on the results
#' pca.crc <- efds %>%
#'   filter_samples(indication == "CRC") %>%
#'   fpcaGadget()
#' report(pca.crc)
#' shine(pca.crc)
#' }
fpcaGadget <- function(x, title = "Principal Components Analysis", ...) {
  frunGadget(fpcaAnalysis, fpcaAnalysisUI, x, title = title, ...)
}

# Embeddable Analysis Module ===================================================

fpcaAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  pca <- callModule(fpcaRun, "pca", rfds, ..., debug = debug)
  view <- callModule(fpcaView, "view", rfds, pca, ..., debug = debug)

  vals <- list(
    result = pca,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacilePcaAnalysisResultContainer",
                   "FacileAnalysisResultContainer")
  vals
}

fpcaAnalysisUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    fpcaRunUI(ns("pca"), debug = debug),
    fpcaViewUI(ns("view"), debug = debug))
}

# Run PCA ======================================================================

#' Minimal shiny module to run fpca
#'
#' @importFrom FacileShine
#'   active_samples
#'   assaySelect
#'   initialized
#'   user
#' @importFrom shiny
#'   callModule
#'   eventReactive
#'   req
fpcaRun <- function(input, output, session, rfds, ..., debug = FALSE) {
  # Provide user with inputs to control:
  # 1. Assay to run PCA on
  # 2. Number of PCs to calculate
  # 3. Number of top varying features to keep
  assert_class(rfds, "ReactiveFacileDataStore")

  assay <- callModule(assaySelect, "assay", rfds)

  result <- eventReactive(input$run, {
    req(initialized(rfds))
    assay_name <- assay$assay_info()$assay
    pcs <- input$pcs
    ntop <- input$ntop
    message("... calculating pca")
    fpca(active_samples(rfds), pcs = pcs, ntop = ntop, assay_name = assay_name,
         custom_key = user(rfds))
  })

  vals <- list(
    result = result,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacilePcaAnalysisResult",
                   "FacilePcaAnalysisResult")
  vals
}

#' @noRd
#' @importFrom FacileShine
#'   assaySelectUI
#' @importFrom shiny
#'   actionButton
#'   column
#'   fluidRow
#'   NS
#'   numericInput
#'   tagList
fpcaRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  # 1. Assay to run PCA on
  # 2. Number of PCs to calculate
  # 3. Number of top varying features to
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), label = "Assay", choices = NULL)),
      column(3, numericInput(ns("pcs"), label = "Number of PCs",
                             value = 10, min = 3, max = 50, step = 1)),
      column(2, numericInput(ns("ntop"), label = "Number of genes",
                             value = 500, min = 50, max = 2000)),
      column(1, actionButton(ns("run"), "Run")))
  )
}

# Visualize PCA ================================================================

#' @noRd
#' @importFrom FacileShine
#'   initialized
#'   categoricalAestheticMap
#' @importFrom plotly
#'   renderPlotly
#' @importFrom shiny
#'   isolate
#'   reactive
#'   req
fpcaView <- function(input, output, session, rfds, result, ..., debug = FALSE) {
  result. <- reactive({
    req(isolate(initialized(rfds)))
    if (is(result, "FacilePCAResult")) {
      out <- result
    } else {
      out <- req(result$result())
    }
    req(is(out, "FacilePCAResult"))
    message("... fpca result retrieved")
    out
  })

  aes <- callModule(categoricalAestheticMap, "aes", rfds,
                    group = FALSE, facet = FALSE, ..., debug = debug)

  output$pcaplot <- renderPlotly({
    res <- req(result.())
    req(is(res, "FacilePCAResult"))

    # NOTE: It looks like covariate() keeps firing, fix!
    # Also FacileShine/tests/shiny-modules/shinytest-facileScatterPlot.R
    # has same issue.
    color. <- aes$color$covariate()
    if (unselected(color.)) color. <- NULL
    shape. <- aes$shape$covariate()
    if (unselected(shape.)) shape. <- NULL

    v <- viz(res, color_aes = color., shape_aes = shape.)
    message("... viz generataed")
    v$plot
  })
}

#' @noRd
#' @importFrom FacileShine categoricalAestheticMapUI
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS
fpcaViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        12,
        wellPanel(categoricalAestheticMapUI(ns("aes"), group = FALSE)))),
    fluidRow(
      column(12, plotlyOutput(ns("pcaplot")))))
}
