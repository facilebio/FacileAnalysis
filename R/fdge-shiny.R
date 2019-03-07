# Full Differential Expression Wrapper =========================================

#' Wrapper module to perform and interact with differential a expression result
#'
#' @export
#' @importFrom shiny callModule
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ShinyFacileDGEResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, with_gsea = FALSE, ...) {
  model <- callModule(fdgeModelDef, "model", rfds)
  dge <- callModule(fdgeRun, "dge", rfds, model, with_gsea = with_gsea, ...)
  interact <- callModule(fdgeViewResult, "view", rfds, dge,
                         with_gsea = with_gsea)

  vals <- list(
    model = model,
    dge = dge,
    interact = interact,
    .ns = session$ns)
  class(vals) <- "ShinyFdgeAnalysis"
  vals
}

#' @noRd
#' @export
#' @importFrom shiny NS tagList tags
fdgeAnalysisUI <- function(id, ...) {
  ns <- NS(id)

  tagList(
    fdgeModelDefUI(ns("model")),
    tags$hr(),
    fdgeRunUI(ns("dge")),
    tags$hr(),
    fdgeViewResultUI(ns("view")))
}

# Just Run fdge ================================================================

#' Shiny module to configure and run fdge over a predefined model.
#'
#' @export
#' @importFrom shiny callModule reactive eventReactive renderText
#'
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
    req(is(out, "FacileDGEModelDefinition"))
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

  class(vals) <- "ShinyFacileDGEResult"
  vals
}

#' @noRd
#' @export
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
fdgeViewResult <- function(input, output, session, rfds, fdge_run, ...) {

  dge.result <- reactive({
    if (is(fdge_run, "FacileDGEResult")) {
      out <- fdge_run
    } else {
      out <- fdge_run$result()
    }
    req(is(out, "FacileDGEResult"))
    out
  })

  dge.stats <- reactive({
    req(dge.result()) %>%
      result() %>%
      select(symbol, feature_id, logFC, FDR = padj, pval) %>%
      arrange(desc(logFC))
  })

  output$stats <- DT::renderDT({
    dat <- req(dge.stats())

    dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
    dtable <- dat %>%
      datatable(filter = "top", extensions = "Scroller", style = "bootstrap",
                class = "compact", width = "100%", rownames = FALSE,
                options = dtopts) %>%
      formatRound(c("logFC", "FDR", "pval"), 3)
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


# Gadget to run a differential expression ======================================

fdgeGadget <- function(dataset, samples = NULL, ...) {
  assert_facile_data_store(dataset)
  if (is.null(samples)) {
    samples <- samples(dataset)
  } else {
    assert_sample_subset(samples, dataset)
  }

  # Invoke the fdgeAnalysis module and UI
}

