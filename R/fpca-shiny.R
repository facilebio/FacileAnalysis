#' Perform a PCA analysis on a (subset) of a FacileDataStore
#'
#' Interactively runs PCA on a FacileDataStore. If you want to run PCA on a
#' subset of samples, pass in a facile_frame sample descriptor
#'
#' @export
#' @examples
#' if (interactive()) {
#' efds <- FacileData::exampleFacileDataSet()
#' # run tumor vs normal comparisons vs each, then run compare9) on the results
#' pca.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
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
#'   batchCorrectConfig
#'   initialized
#'   user
#' @importFrom shiny
#'   callModule
#'   eventReactive
#'   req
#'   updateNumericInput
#'   withProgress
#' @importFrom shinyjs toggleState
fpcaRun <- function(input, output, session, rfds, ..., debug = FALSE,
                    .reactive = TRUE) {

  # Provide user with inputs to control:
  # 1. Assay to run PCA on
  # 2. Number of PCs to calculate
  # 3. Number of top varying features to keep
  assert_class(rfds, "ReactiveFacileDataStore")

  isolate. <- if (.reactive) base::identity else shiny::isolate

  active.samples <- reactive({
    req(initialized(rfds))
    # isolate.(active_samples(rfds))
    ftrace("Updating active samples")
    active_samples(rfds)
  })

  batch <- callModule(batchCorrectConfig, "batch", rfds)

  # Set the maximum number of PCs such that they do not exceed number of samples
  # Running PCA is also disabled if there are less than three samples.
  observeEvent(active.samples(), {
    asamples. <- req(active.samples())
    nsamples <- nrow(asamples.)
    value <- if (input$pcs > nsamples) nsamples else NULL
    updateNumericInput(session, "pcs", value = value, max = nsamples)
    toggleState("run", condition = nsamples >= 3L)
  })

  assay <- callModule(assaySelect, "assay", rfds)

  result <- eventReactive(input$run, {
    req(initialized(rfds))
    samples. <- active.samples()
    assay_name <- assay$assay_info()$assay
    pcs <- input$pcs
    ntop <- input$ntop
    batch. <- name(batch$batch)
    main. <- name(batch$main)
    withProgress({
      fpca(samples., dims = pcs, ntop = ntop, assay_name = assay_name,
           batch = batch., main = main., custom_key = user(rfds))
    }, message = "Performing PCA")
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
#'   batchCorrectConfigUI
#' @importFrom shiny
#'   actionButton
#'   column
#'   fluidRow
#'   NS
#'   numericInput
#'   tagList
#' @importFrom shinyWidgets dropdownButton
fpcaRunUI <- function(id, width_opts = "200px", ..., debug = FALSE) {
  ns <- NS(id)
  # 1. Assay to run PCA on
  # 2. Number of PCs to calculate
  # 3. Number of top varying features to
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), label = "Assay", choices = NULL)),
      column(4, batchCorrectConfigUI(ns("batch"), direction = "horizontal")),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          dropdownButton(
            inputId = ns("opts"),
            icon = icon("sliders"),
            status = "primary", circle = FALSE,
            width = width_opts,

            numericInput(ns("pcs"), label = "Number of PCs",
                         value = 10, min = 2, max = 30, step = 1),
            numericInput(ns("ntop"), label = "Number of genes",
                         value = 500, min = 50, max = 5000, step = 500)))),
      column(1, actionButton(ns("run"), "Run"), style = "margin-top: 1.7em")
    ))
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
    color_aes = NULL,
    shape_aes = NULL,
    hover = NULL,
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
    pcs <- setNames(seq(pcs),
                    sprintf("PC%d (%0.1f%% variance)", as.integer(pcs),
                            pca.$percent_var * 100))

    updateSelectInput(session, "xaxis", choices = pcs, selected = pcs[1])
    updateSelectInput(session, "yaxis", choices = pcs, selected = pcs[2])
    updateSelectInput(session, "zaxis", choices = c("---", pcs), selected = "")

    pc.choices <- setNames(
      names(pca.$percent_var),
      sprintf("%s (%0.1f%%)", names(pca.$percent_var), pca.$percent_var * 100))

    updateSelectInput(session, "loadingsPC", choices = pcs, selected = pcs[1L])
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

  # Color Mapping ..............................................................
  # This will be painful. Let's enable the PCA plot to be colored by a
  # categorical covariate, or one of the features picked from the loadings on
  # the feature table

  # Color covariate was selected from the categoricalAestheticMap module
  observeEvent(aes$map(), {
    aes.map <- aes$map()
    acolor <- aes.map$color
    ashape <- aes.map$shape
    ahover <- aes.map$hover
    if (!identical(state$color_aes, acolor)) state$color_aes <- acolor
    if (!identical(state$shape_aes, ashape)) state$shape_aes <- ashape
    if (!identical(state$hover, ahover)) state$hover <- ahover
  })

  # A gene was selected from the PC loadings table
  observeEvent(input$loadings_rows_selected, {
    selected <- input$loadings_rows_selected
    if (!is.null(selected)) {
      dat <- pc.loadings()
      feature <- paste0("feature:", dat[["feature_id"]][selected])
      if (!identical(state$color_aes, feature)) {
        state$color_aes <- feature
        # TODO: update selected color selection in categoricalAestheticMap
        # `aes` module
      }
    }
  })

  # Render PCA Plot ............................................................
  pcaviz <- reactive({
    pca. <- req(pca())
    pc.calcd <- pcs_calculated()
    req(length(pc.calcd) >= 2)
    axes <- intersect(c(input$xaxis, input$yaxis, input$zaxis), seq(pc.calcd))
    req(length(axes) >= 2)

    aes.map <- aes$map()
    acolor <- state$color_aes
    ashape <- state$shape_aes
    ahover <- state$hover

    if (identical(substr(acolor, 1L, 8L), "feature:")) {
      feature_id <- sub("feature:", "", state$color_aes)
      current.cols <- colnames(pca.$result)
      pca.$result <- with_assay_data(pca.$result, features = feature_id,
                                     normalized = TRUE,
                                     batch = param(pca., "batch"),
                                     main = param(pca., "main"))
      added <- setdiff(colnames(pca.$result), current.cols)
      acolor <- added
    }

    viz(pca., axes, color_aes = acolor, shape_aes = ashape, hover = ahover,
        width = NULL, height = 550)
  })

  output$pcaplot <- renderPlotly({
    plot(req(pcaviz()))
  })

  feature.ranks <- reactive({
    req(pca()) %>%
      signature(signed = TRUE, ntop = 50) %>%
      tidy()
  })

  pc.loadings <- reactive({
    pc <- paste0("PC", input$loadingsPC)
    req(feature.ranks()) %>%
      filter(dimension == pc) %>%
      select(symbol, feature_id, score)
  })

  # Loadings Table .............................................................
  output$loadings <- DT::renderDT({
    dtopts <- list(deferRender = TRUE, scrollY = 450,
                   # scroller = TRUE,
                   pageLength = 15,
                   lengthMenu = c(15, 30, 50))

    pc.dat <- pc.loadings()
    num.cols <- colnames(pc.dat)[sapply(pc.dat, is.numeric)]

    dt <- datatable(pc.dat, filter = "top",
                    style = "bootstrap",
                    class = "display", width = "100%", rownames = FALSE,
                    selection = "single",
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
    tags$h4(
      tags$span("PCA Result", style = "background: #fff; padding: 0 10px 0 0"),
      style = "border-bottom: 1px solid #000; line-height: 0.1em; margin-bottom: 13px"),
    fluidRow(
      column(
        7,
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
        plotlyOutput(ns("pcaplot"))),
      column(
        5,
        tags$h4("Feature Loadings"),
        tags$div(selectInput(ns("loadingsPC"), "Principal Component",
                             choices = NULL)),
        withSpinner(DT::DTOutput(ns("loadings"), height = "650px")))))
}
