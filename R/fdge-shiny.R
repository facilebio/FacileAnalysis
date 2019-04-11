# Full Differential Expression Wrapper =========================================

#' Wrapper module to perform and interact with differential a expression result
#'
#' @export
#' @importFrom shiny callModule
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ShinyFacileDGEResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, ...) {
  model <- callModule(fdgeModelDef, "model", rfds, ...)
  dge <- callModule(fdgeRun, "dge", rfds, model, ...)
  view <- callModule(fdgeViewResult, "view", rdfs, dge, ...)

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
#' @importFrom shiny NS tagList tags
#' @importFrom multiGSEA failWith
fdgeAnalysisUI <- function(id, ...) {
  ns <- NS(id)

  tagList(
    fdgeModelDefUI(ns("model")),
    tags$hr(),
    fdgeRunUI(ns("dge")),
    tags$hr(),
    fdgeViewResultUI(ns("view")))
}

#' @noRd
#' @export
#' @importFrom shiny browserViewer dialogViewer
#' @examples
#' if (ineteractive()) {
#'   FacileData::exampleFacileDataSet() %>%
#'     filter_samples(indication == "CRC") %>%
#'     fdgeGadget()
#' }
fdgeGadget <- function(x, user = Sys.getenv("USER"),
                       title = "Differential Expression Analysis",
                       height = 600, width = 800, ..., debug = FALSE) {
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
    miniContentPanel(fdgeAnalysisUI("analysis")),
    NULL)

  server <- function(input, output, session) {
    rfds <- callModule(reactiveFacileDataStore, "ds", fds., samples., user)
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)

    observeEvent(input$done, {
      result. <- failWith(NULL, analysis$fdge$result())
      annotation <- tibble(
        # TODO: Get annotations from brushing on volcano?
        feature_type = character(),
        feature_id = character(),
        annotation = character())
      out <- list(
        result = result.,
        annotation = annotation)
      class(out) <- c("FacileDGEGadgetResult", "FacileGadgetResult",
                      "FacileAnalysisResult")
      stopApp(invisible(out))
    })
    observeEvent(input$cancel, {
      stopApp(invisible(NULL))
    })
  }

  if (debug) {
    viewer <- browserViewer()
  } else {
    viewer <- dialogViewer(title, height = height, width = width)
  }

  runGadget(ui, server, viewer = viewer, stopOnCancel = FALSE)
}

# fgadget <- function(x, module, ...) {
#   UseMethod("fgadget", x)
# }


# Just Run fdge ================================================================

#' Shiny module to configure and run fdge over a predefined model.
#'
#' @export
#' @importFrom shiny callModule reactive eventReactive renderText
#' @importFrom FacileShine
#'   assaySelect
#'   unselected
#' @param model A linearm model definition. Can be either a "naked"
#'   `FacileDGEModelDefinition` that is returned from a call to
#'   [fdge_model_def()], or the `ShinyDGEModelDefinition` object returned from
#'   the [fdgeModelDef()] module.
#' @param with_gsea Include option to run a GSEA?
#' @param ... passed into [fdge()]
#' @return A list of stuff. `$result` holds the `FacileDGEResult` wrapped in
#'   a `reactive()`, ie. a `ReactiveDGEResult`.
fdgeRun <- function(input, output, session, rfds, model, with_gsea = FALSE, ...,
                    .reactive = TRUE) {
  kosher <- is(model, "FacileDGEModelDefinition") ||
    is(model, "ShinyDGEModelDefinition")
  if (!kosher) {
    stop("Invalid object passed as `model`: ", class(model)[1L])
  }

  assay <- callModule(assaySelect, "assay", rfds, .reactive = .reactive)

  rmodel <- reactive({
    out <- if (is(model, "FacileDGEModelDefinition")) model else model$result()
    req(is(out, "FacileDGEModelDefinition"),
        !is(out, "IncompleteModelDefintion"))
    out
  })

  dge <- eventReactive(input$run, {
    model. <- req(rmodel())
    aname <- assay$assay_info()$assay
    req(!unselected(aname))
    fdge(model., assay_name = aname)
  })

  output$debug <- shiny::renderText({
    dge. <- req(dge())
    format(dge.)
  })

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
#'   column
#'   fluidRow
#'   NS
#'   tagList
#'   tags
#'   verbatimTextOutput
fdgeRunUI <- function(id, ...) {
  ns <- NS(id)
  # Need widgets for:
  # 1. assay_name selection
  # 2. method selection (voom, qlf, etc.)
  # 3. treat_lfc
  # 4. filter criterion (let's leave this blank for now.)
  # 5. Run button

  out <- tagList(
    fluidRow(
      column(1, assaySelectUI(ns("assay"), label = "Assay", choices = NULL)),
      column(1, actionButton(ns("run"), "Run"))
    ),
    NULL)

  out <- tagList(
    out,
    tags$hr(),
    shiny::verbatimTextOutput(ns("debug"), placeholder = TRUE)
  )
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
fdgeViewResult <- function(input, output, session, rfds, fdge_result, ...) {

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
fdgeViewResultUI <- function(id, ...) {
  ns <- NS(id)
  DT::DTOutput(ns("stats"))
}
