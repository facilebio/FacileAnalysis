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
fpcaGadget <- function(x, title = "Principal Components Analysis",
                       viewer = "browser", ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fpcaAnalysis, fpcaAnalysisUI, x, title = title, viewer = viewer,
             ...)
}

# Embeddable Analysis Module ===================================================

#' @noRd
#' @export
fpcaAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  pca <- callModule(fpcaRun, "pca", rfds, ..., debug = debug)
  view <- callModule(fpcaView, "view", rfds, pca, ..., debug = debug)

  # Only show view widgets when there is a pca result
  observe({
    res. <- req(faro(pca))
    show <- is(res., "FacilePcaAnalysisResult")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    main = pca,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacilePcaAnalysisResultContainer",
                   "ReactiveFacileAnalysisResultContainer")
  vals
}

#' @noRd
#' @export
#' @importFrom shinyjs hidden
fpcaAnalysisUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    tags$div(id = ns("runbox"), fpcaRunUI(ns("pca"), debug = debug)),
    hidden(tags$div(id = ns("viewbox"), fpcaViewUI(ns("view"), debug = debug))))
}

# Run PCA ======================================================================

#' Minimal shiny module to run fpca
#'
#' @export
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
    fpca(active_samples(rfds), dims = pcs, ntop = ntop, assay_name = assay_name,
         custom_key = user(rfds))
  })

  vals <- list(
    faro = result,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacilePcaAnalysisResult",
                   "ReactiveFacileAnalysisResult",
                   "FacilePcaAnalysisResult")
  vals
}

#' @noRd
#' @export
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
#' @export
#' @importFrom FacileShine
#'   initialized
#'   categoricalAestheticMap
#' @importFrom FacileViz with_aesthetics fscatterplot input_data
#' @importFrom plotly
#'   renderPlotly
#' @importFrom shiny
#'   isolate
#'   reactive
#'   req
fpcaView <- function(input, output, session, rfds, pcares, ...,
                     feature_selection = session$ns("features"),
                     sample_selection = session$ns("samples"),
                     debug = FALSE) {

  state <- reactiveValues(
    # store brushed samples from scatterplot
    scatter_select = tibble(assay_name = character(), feature_id = character()))

  pca <- reactive({
    req(initialized(pcares))
    faro(pcares)
  })

  pcs_calculated <- reactive({
    pca. <- req(pca())
    names(pca.$percent_var)
  })

  # When a new fpca result is produced, we may want to reset a few things.
  # Trigger that UI restting/cleanup work here.
  observeEvent(pca(), {
    pca. <- req(pca())
    pcs <- req(pcs_calculated())
    pcs <- setNames(seq(pcs), pcs)

    updateSelectInput(session, "xaxis", choices = pcs, selected = pcs[1])
    updateSelectInput(session, "yaxis", choices = pcs, selected = pcs[2])
    updateSelectInput(session, "zaxis", choices = pcs, selected = "")
  }, priority = 5) # upping priority so some withProgress things hide quick

  assay_name. <- reactive(param(req(pca()), "assay_name"))
  samples. <- reactive(samples(req(pca())))

  batch_corrected <- reactive({
    batch <- param(req(pca()), "batch")
    !is.null(batch)
  })

  aes <- callModule(categoricalAestheticMap, "aes", rfds,
                    color = TRUE, shape = TRUE, group = FALSE, facet = FALSE,
                    hover = TRUE, ..., debug = debug)

  pcaviz <- reactive({
    pca. <- req(pca())
    pc.calcd <- pcs_calculated()
    req(length(pc.calcd) >= 2)
    axes <- intersect(c(input$xaxis, input$yaxis, input$zaxis), seq(pc.calcd))
    req(length(axes) >= 2)

    aes.map <- aes$map()
    viz(pca., axes, color_aes = aes.map$color, shape_aes = aes.map$shape,
        hover = aes.map$hover, width = NULL, height = NULL)
  })

  output$pcaplot <- renderPlotly({
    plot(req(pcaviz()))
  })

  feature.ranks <- reactive({
    req(pca()) %>%
      signature(signed = TRUE, ntop = 50) %>%
      tidy()
  })

  output$loadings <- DT::renderDT({
    dat <- req(feature.ranks()) %>%
      select(dimension, symbol, feature_id, score) %>%
      mutate(dimension = factor(dimension))
    dtopts <- list(deferRender = TRUE, scrollY = 300,
                   # scroller = TRUE,
                   pageLength = 15,
                   lengthMenu = c(15, 30, 50))
    num.cols <- colnames(dat)[sapply(dat, is.numeric)]

    dt <- datatable(dat, filter = "top",
                    style = "bootstrap",
                    class = "display", width = "100%", rownames = FALSE,
                    selection = "none",
                    options = dtopts)
    formatRound(dt, num.cols, 3)
  }, server = TRUE)
}

#' @noRd
#' @export
#' @importFrom FacileShine categoricalAestheticMapUI
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS
#' @importFrom shinycssloaders withSpinner
fpcaViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, selectInput(ns("xaxis"), "X axis", choices = NULL)),
      column(4, selectInput(ns("yaxis"), "Y axis", choices = NULL)),
      column(4, selectInput(ns("zaxis"), "Z axis", choices = NULL))),
    fluidRow(
      column(
        12,
        wellPanel(
          categoricalAestheticMapUI(
            ns("aes"),
            color = TRUE, shape = TRUE, hover = TRUE,
            group = FALSE)))),
    fluidRow(
      column(7, plotlyOutput(ns("pcaplot"))),
      column(
        5,
        tags$h4("Feature Loadings"),
        withSpinner(DT::DTOutput(ns("loadings"))))))
}
