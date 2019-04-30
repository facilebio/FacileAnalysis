
#' @noRd
#' @export
#' @importFrom FacileShine ReactiveFacileDataStore
#' @importFrom shiny
#'   browserViewer
#'   dialogViewer
#' @importFrom multiGSEA
#'   failWith
#' @examples
#' if (interactive()) {
#' efds <- exampleFacileDataSet()
#' # run tumor vs normal comparisons vs each, then run compare9) on the results
#' dge.crc <- efds %>%
#'   filter_samples(indication == "CRC") %>%
#'   fdgeGadget()
#' dge.blca <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdgeGadget()
#' dge.comp <- compare(dge.crc, dge.blca)
#' if (interactive()) {
#'   report(dge.comp)
#'   shine(dge.comp)
#' }
fdgeGadget <- function(x, user = Sys.getenv("USER"),
                       title = "Differential Expression Analysis",
                       height = 600, width = 1000,
                       viewer = "dialog", ...,
                       # viewer = "browser", ...,
                       debug = FALSE) {
  bs4dash <- getOption("facile.bs4dash")
  options(facile.bs4dash = FALSE)
  on.exit(options(facile.bs4dash = bs4dash))

  viewer <- gadget_viewer(viewer, title, width, height)

  if (is(x, "facile_frame")) {
    fds. <- fds(x)
    samples. <- x
  } else {
    fds. <- x
    samples. <- samples(x)
  }

  assert_class(fds., "FacileDataStore")
  assert_class(samples., "facile_frame")

  ui <- miniPage(
    gadgetTitleBar(title),
    miniContentPanel(fdgeAnalysisUI("analysis", debug = debug)),
    NULL)

  server <- function(input, output, session) {
    # rfds <- callModule(reactiveFacileDataStore, "ds", fds., samples., user)
    rfds <- ReactiveFacileDataStore(fds., "ds", user = user, samples = samples.)
    analysis <- callModule(fdgeAnalysis, "analysis", rfds, debug = debug)

    observeEvent(input$done, {
      result. <- failWith(NULL, analysis$fdge$result())
      result.[["fds"]] <- fds.
      result.[["params"]][["model_def"]][["fds"]] <- fds.

      annotation <- FacileShine:::.empty_feature_annotation_tbl()
      out <- list(
        result = result.,
        annotation = annotation,
        fds = fds.)
      class(out) <- c("FacileDGEGadgetResult",
                      "FacileGadgetResult",
                      "FacileAnalysisResult")
      stopApp(invisible(out))
    })
    observeEvent(input$cancel, {
      stopApp(invisible(NULL))
    })
  }

  runGadget(ui, server, viewer = viewer, stopOnCancel = FALSE)
}

#' @noRd
report.FacileDGEGadgetResult <- function(x, ...) {
  report(result(x), ...)
}

#' @noRd
viz.FacileDGEGadgetResult <- function(x, ...) {
  viz(result(x), ...)
}

#' @noRd
shine.FacileDGEGadgetResult <- function(x, ...) {
  shine(result(x), ...)
}

#' @noRd
compare.FacileDGEGadgetResult <- function(x, y, ...) {
  if (is(y, "FacileDGEGadgetResult")) y <- result(y)
  compare(result(x), y)
}

#' Wrapper module to perform and interact with a differential expression result.
#'
#' This module can be embedded within a shiny app, or called from a gadget.
#'
#' @export
#' @importFrom shiny callModule
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ShinyFacileDGEResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  model <- callModule(fdgeModelDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)
  view <- callModule(fdgeView, "view", rdfs, dge, ..., debug = debug)

  vals <- list(
    fdge = dge,
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
#' @importFrom bs4Dash
#'   bs4Box
#' @importFrom shinydashboard
#'   box
fdgeAnalysisUI <- function(id, ..., debug = FALSE,
                           bs4dash = isTRUE(getOption("facile.bs4dash"))) {
  ns <- NS(id)
  box. <- if (bs4dash) bs4Dash::bs4Box else shinydashboard::box
  tagList(
    box.(title = "Model Definition", width = 12,
         fdgeModelDefRunUI(ns("model"), debug = debug)),
    box.(title = "Testing Parameters", width = 12,
         fdgeRunUI(ns("dge"), debug = debug)),
    box.(title = "Testing Results", width = 12,
         fdgeViewUI(ns("view"), debug = debug)))
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
    req(is(out, "FacileDGEModelDefinition"),
        !is(out, "IncompleteModelDefintion"))
    out
  })

  assay <- callModule(assaySelect, "assay", rfds, .reactive = .reactive)

  # Update the dge_methods available given the selected assay
  observe({
    assay_type. <- assay$assay_info()$assay_type
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

  dge <- eventReactive(input$run, {
    model. <- req(rmodel())
    assay_name. <- assay$assay_info()$assay
    method. <- input$dge_method
    req(
      !unselected(assay_name.),
      !unselected(method.))
    sample_weights. <- sample_weights()
    treat_lfc. <- log2(input$treatfc)
    fdge(model., assay_name = assay_name., method = method.,
         with_sample_weights = sample_weights., treat_lfc = treat_lfc.)
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
#'   NS
#'   numericInput
#'   selectInput
#'   tagList
#'   tags
#'   verbatimTextOutput
fdgeRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), label = "Assay", choices = NULL)),
      column(3, selectInput(ns("dge_method"), label = "Method", choices = NULL)),
      column(3, numericInput(ns("treatfc"), label = "Min. Fold Change",
                             value = 1, min = 1, max = 10, step = 0.25)),
      column(2, checkboxInput(ns("weights"), label = "Sample Weights")),
      column(1, actionButton(ns("run"), "Run")))
    )

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
#' @importFrom DT datatable renderDT formatRound
#' @importFrom shiny req
#'
#' @param fdge_result The result from calling [fdge()].
fdgeView <- function(input, output, session, rfds, fdge_result, ...,
                     debug = FALSE) {

  # The FacileDGEResult object
  result. <- reactive({
    if (is(fdge_result, "FacileDGEResult")) {
      out <- fdge_result
    } else {
      out <- fdge_result$result()
    }
    req(is(out, "FacileDGEResult"))
    out
  })

  dge.stats <- reactive({
    .fdge <- req(result.())
    .result <- ranks(.fdge) %>% result()
    if (.fdge[["test_type"]] == "ttest") {
      out <- select(.result, symbol, feature_id, logFC, FDR = padj, pval)
    } else {
      out <- select(.result, symbol, feature_id, F, FDR = padj, pval)
    }
    out
  })

  output$stats <- DT::renderDT({
    dat <- req(dge.stats())
    num.cols <- colnames(dat)[sapply(dat, is.numeric)]

    dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
    dtable <- dat %>%
      datatable(filter = "top", extensions = "Scroller", style = "bootstrap",
                class = "compact", width = "100%", rownames = FALSE,
                options = dtopts) %>%
      formatRound(num.cols, 3)
    dtable
  }, server = TRUE)
}

#' @noRd
#' @export
#' @importFrom shiny NS
#' @importFrom DT DTOutput
fdgeViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  DT::DTOutput(ns("stats"))
}
