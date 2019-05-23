
#' @noRd
#' @export
#' @examples
#' if (interactive()) {
#' # run tumor vs normal comparisons vs each, then run compare() on the results
#' options(facile.log.level.fshine = "trace")
#' efds <- exampleFacileDataSet()
#' dge.crc <- efds %>%
#'   filter_samples(indication == "CRC") %>%
#'   fdgeGadget(viewer = "pane")
#' dge.blca <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdgeGadget(viewer = "pane")
#' dge.comp <- compare(dge.crc, dge.blca)
#'
#' report(dge.comp)
#' shine(dge.comp)
#' }
fdgeGadget <- function(x, title = "Differential Expression Analysis",
                       height = 800, width = 1000, ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fdgeAnalysis, fdgeAnalysisUI, x, title = title,
             height = height, width = width, ...)
}

#' Wrapper module to perform and interact with a differential expression result.
#'
#' This module can be embedded within a shiny app, or called from a gadget.
#'
#' @export
#' @importFrom shiny callModule
#' @importFrom shinyjs toggleElement
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ShinyFacileDGEResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  model <- callModule(fdgeModelDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)
  view <- callModule(fdgeView, "view", rfds, dge, ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = debug)

  # Only show the view UI when there is (at least) a defined model. When
  # done this way, the showSpinner() around the view's output datatable works
  # like a "wait, processing DGE" while it's running.
  observe({
    model. <- model$result()
    show <- is(model., "FacileDGEModelDefinition")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    result = dge,
    model = model,
    view = view,
    .ns = session$ns)
  class(vals) <- "ShinyFdgeAnalysis"
  vals
}

#' @noRd
#' @export
#' @importFrom shiny
#'   fluidRow
#'   NS
#'   tagList
#'   tags
#'   wellPanel
#' @importFrom shinydashboard box
#' @importFrom shinyjs hidden
fdgeAnalysisUI <- function(id, ..., debug = FALSE,
                           bs4dash = isTRUE(getOption("facile.bs4dash"))) {
  ns <- NS(id)
  # box. <- if (bs4dash) bs4Dash::bs4Box else shinydashboard::box
  box. <- shinydashboard::box
  tagList(
    tags$div(id = ns("modelbox"),
             box.(title = "Model Definition", width = 12,
                  fdgeModelDefRunUI(ns("model"), debug = debug))),
    tags$div(id = ns("dgebox"),
             box.(title = "Testing Parameters", width = 12,
                  fdgeRunUI(ns("dge"), debug = debug))),
    hidden(tags$div(id = ns("viewbox"),
                    box.(title = "Testing Results", width = 12,
                         fdgeViewUI(ns("view"), debug = debug)))))
}

# Building block fdge shiny modules ============================================

#' Shiny module to configure and run fdge over a predefined model.
#'
#' @export
#' @importFrom shiny
#'   callModule
#'   eventReactive
#'   reactive
#'   renderText
#'   withProgress
#' @importFrom shinyjs toggleState
#' @importFrom FacileShine
#'   assaySelect
#'   unselected
#' @param model A linear model definition. Can be either a "naked"
#'   `FacileDGEModelDefinition` that is returned from a call to
#'   [fdge_model_def()], or the `ShinyDGEModelDefinition` object returned from
#'   the [fdgeModelDefRun()] module.
#' @param with_gsea Include option to run a GSEA?
#' @param ... passed into [fdge()]
#' @return A list of stuff. `$result` holds the `FacileDGEResult` wrapped in
#'   a `reactive()`, ie. a `ReactiveDGEResult`.
fdgeRun <- function(input, output, session, rfds, model, with_gsea = FALSE, ...,
                    debug = FALSE, .reactive = TRUE) {
  kosher <- is(model, "FacileDGEModelDefinition") ||
    is(model, "ShinyDGEModelDefinition")
  if (!kosher) {
    stop("Invalid object passed as `model`: ", class(model)[1L])
  }

  rmodel <- reactive({
    out <- if (is(model, "FacileDGEModelDefinition")) model else model$result()
    # req(is(out, "FacileDGEModelDefinition"),
    #     !is(out, "IncompleteModelDefintion"))
    out
  })

  assay <- callModule(assaySelect, "assay", rfds, .reactive = .reactive)

  # Update the dge_methods available given the selected assay
  observe({
    ainfo <- req(assay$assay_info())
    assay_type. <- ainfo$assay_type


    req(!unselected(assay_type.))
    method. <- input$dge_method
    methods. <- fdge_methods(assay_type.)$dge_method
    selected. <- if (method. %in% methods.) method. else methods.[1L]
    updateSelectInput(session, "dge_method", choices = methods.,
                      selected = selected.)
  })

  sample_weights <- reactive({
    method. <- input$dge_method
    req(!unselected(method.))
    can_weight <- fdge_methods() %>%
      filter(dge_method == method.) %>%
      pull(can_sample_weight)
    input$weights && can_weight
  })

  observe({
    model. <- rmodel()
    enable <- is(model., "FacileDGEModelDefinition") &&
      !is(model., "FacileFailedModelDefinition") &&
      !is(model., "IncompleteModelDefintion")
    toggleState("run", condition = enable)
  })

  dge <- eventReactive(input$run, {
    model. <- rmodel()
    cando <- is(model., "FacileDGEModelDefinition") &&
      !is(model., "FacileFailedModelDefinition") &&
      !is(model., "IncompleteModelDefintion")
    req(cando)
    assay_name. <- assay$assay_info()$assay
    method. <- input$dge_method
    req(
      !unselected(assay_name.),
      !unselected(method.))
    sample_weights. <- sample_weights()
    treat_lfc. <- log2(input$treatfc)

    withProgress({
      fdge(model., assay_name = assay_name., method = method.,
           with_sample_weights = sample_weights., treat_lfc = treat_lfc.)
    }, message = "Performing differential expression")
  })

  if (debug) {
    output$debug <- shiny::renderText({
      dge. <- req(dge())
      format(dge.)
    })
  }

  vals <- list(
    result = dge,
    .ns = session$ns)

  class(vals) <- "FacileShinyDGEResult"
  vals
}

#' @noRd
#' @export
#' @importFrom FacileShine
#'   assaySelectUI
#' @importFrom shiny
#'   actionButton
#'   checkboxInput
#'   column
#'   fluidRow
#'   icon
#'   NS
#'   numericInput
#'   selectInput
#'   tagList
#'   tags
#'   verbatimTextOutput
#' @importFrom shinyWidgets dropdownButton prettyCheckbox
fdgeRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), choices = NULL)),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          dropdownButton(
            inputId = ns("opts"),
            icon = icon("sliders"),
            status = "primary", circle = FALSE,
            width = "300px",
            selectInput(ns("dge_method"), label = "Method", choices = NULL),
            numericInput(ns("treatfc"), label = "Min. Fold Change",
                         value = 1, min = 1, max = 10, step = 0.25),
            prettyCheckbox(ns("weights"), label = "Use Sample Weights",
                           status = "primary")))),
      column(1, actionButton(ns("run"), "Run", style = "margin-top: 1.7em"))
    ))

  if (debug) {
    out <- tagList(
      out,
      tags$hr(),
      verbatimTextOutput(ns("debug"), placeholder = TRUE))
  }

  out
}

# Visualize the fdge result ====================================================

#' Shiny module that interacts with result of fdgeRun module
#'
#' Provides an interactive view over the result of running the [fdgeRun()]
#' module.
#'
#' @export
#' @importFrom FacileViz fscatterplot
#' @importFrom DT datatable dataTableProxy formatRound renderDT replaceData
#' @importFrom plotly event_data
#' @importFrom shiny req
#' @param result The result from calling [fdge()].
fdgeView <- function(input, output, session, rfds, result, ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = FALSE) {

  # The FacileDGEResult object
  result. <- reactive({
    req(isolate(initialized(rfds)))
    if (is(result, "FacileDGEResult")) {
      out <- result
    } else {
      out <- result$result()
    }
    req(is(out, "FacileDGEResult"))
    out
  })

  dge.stats.all <- reactive({
    .fdge <- req(result.())
    .result <- ranks(.fdge) %>% result()
    if (.fdge[["test_type"]] == "ttest") {
      out <- select(.result, symbol, feature_id, logFC, FDR = padj, pval)
    } else {
      out <- select(.result, symbol, feature_id, F, FDR = padj, pval)
    }
    out
  })

  output$volcano <- renderPlotly({
    dat. <- req(dge.stats.all())
    dat.$yaxis <- -log10(dat.$pval)
    fplot <- fscatterplot(dat., c("logFC", "yaxis"),
                          xlabel = "logFC", ylabel = "-log10(pval)",
                          width = NULL, height = NULL,
                          hover = c("symbol", "logFC", "FDR"),
                          webgl = TRUE, event_source = feature_selection,
                          key = "feature_id")
    fplot$plot
  })

  # NULL or a character vector of feature_ids that are selected
  selected_features <- reactive({
    selected <- event_data("plotly_selected", source = feature_selection)
    if (!is.null(selected)) selected <- selected$key
    selected
  })

  dge.stats <- reactive({
    stats. <- req(dge.stats.all())
    selected <- selected_features()
    if (!is.null(selected)) stats. <- filter(stats., feature_id %in% selected)
    stats.
  })

  output$boxplot <- renderPlotly({

  })

  selected_samples <- reactive({
    # return the brushed samples from boxplot
    NULL
  })

  output$stats <- DT::renderDT({
    # Triggered when new result comes in
    res <- req(result.())
    dat <- req(dge.stats())
    num.cols <- colnames(dat)[sapply(dat, is.numeric)]

    dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
    dtopts <- list()

    dtable <- dat %>%
      datatable(filter = "top",
                # extensions = "Scroller",
                style = "bootstrap",
                class = "display", width = "100%", rownames = FALSE,
                selection = "single",
                options = dtopts) %>%
      formatRound(num.cols, 3)
    dtable
  }, server = TRUE)
  # dtproxy <- dataTableProxy("stats")
  # observe({
  #   stats. <- req(dge.stats())
  #   replaceData(dtproxy, stats., resetPaging = FALSE)
  # })

  vals <- list(
    selected_features = selected_features,
    selected_samples = selected_samples)
  return(vals)
}

#' @noRd
#' @export
#' @importFrom DT DTOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shiny column fluidRow NS
#' @importFrom shinycssloaders withSpinner
fdgeViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  fluidRow(
    column(
      4,
      # plotlyOutput(ns("boxplot")),
      withSpinner(plotlyOutput(ns("volcano")))),
    column(8, withSpinner(DT::DTOutput(ns("stats"))))
  )

}
