#' @include fdge.R
NULL

# shiny modules to tweak the options required to run an fdge analysis

# Method options ===============================================================

#' @noRd
#' @importFrom shinyWidgets updatePrettyCheckbox updatePickerInput
#' @param rfds ReactiveFacileDataStore
#' @param model a FacileLinearModelDefinition (or ReactiveFacileLinearModelDefinition)
fdgeRunOptions <- function(input, output, session, rfds, model, assay, ...) {
  kosher <- is(model, "FacileLinearModelDefinition") ||
    is(model, "ReactiveFacileLinearModelDefinition")
  if (!kosher) {
    stop("Invalid object passed as `model`: ", class(model)[1L])
  }

  state <- reactiveValues(
    dge_method = "__initializing__",
    feature_class = "__initializing__")

  # When the assay type changes, we need to update the type of testing we can do.
  observeEvent(assay$assay_info(), {
    ainfo <- req(assay$assay_info())
    assay_type. <- req(ainfo$assay_type, !unselected(ainfo$assay_type))

    methods. <- fdge_methods(assay_type., on_missing = "warning")$dge_method
    selected <- state$dge_method
    if (!selected %in% methods.) selected <- methods.[1L]
    updateSelectInput(session, "dge_method", choices = methods.,
                      selected = selected)

    features. <- assay$features()
    ftypes <- unique(features.[["meta"]])
    if (ainfo$assay_type == "rnaseq" && "protein_coding" %in% ftypes) {
      selected <- "protein_coding"
    } else {
      selected <- ftypes
    }
    updatePickerInput(session, "features", choices = ftypes,
                      selected = selected)
  })

  # Sets up appropriate UI to retreive the params used to filter out
  # features measured using the selected assay
  # ffilter <- callModule(fdgeFeatureFilter, "ffilter", rfds, model,
  #                       assay, ..., debug = debug)
  ffilter <- list(
    universe = reactive({
      filter(assay$features(), .data$meta %in% input$features)
    })
  )

  observeEvent(input$dge_method, {
    if (input$dge_method != state$dge_method) {
      state$dge_method <- input$dge_method
    }
  })

  dge_method <- reactive(state$dge_method)

  can_sample_weight <- reactive({
    method. <- req(dge_method(), !unselected(dge_method()))
    ainfo <- req(assay$assay_info())
    atype <- req(ainfo$assay_type, !unselected(ainfo$assay_type))
    fdge_methods() %>%
      filter(.data$assay_type == .env$atype, .data$dge_method == .env$method.) %>%
      pull(can_sample_weight) %>%
      isTRUE()
  })

  observeEvent(can_sample_weight(), {
    can.weight <- can_sample_weight()
    if (!can.weight) {
      updatePrettyCheckbox(session, "sample_weights", value = FALSE)
    }
    toggleState("sample_weights", condition = can.weight)
  })

  sample_weights <- reactive({
    isTRUE(can_sample_weight() && input$sample_weights)
  })

  treat_lfc <- reactive(log2(input$treatfc))

  vals <- list(
    dge_method = dge_method,
    sample_weights = sample_weights,
    treat_lfc = treat_lfc,
    feature_filter = ffilter,
    .state = state,
    .ns = session$ns)

  class(vals) <- c("FdgeRunOptions")
  vals
}

#' @noRd
#' @importFrom shinyWidgets dropdown prettyCheckbox pickerInput
#' @importFrom shiny icon wellPanel
fdgeRunOptionsUI <- function(id, width = "300px", ..., debug = FALSE) {
  ns <- NS(id)
  dropdown(
    inputId = ns("opts"),
    icon = icon("sliders"),
    status = "primary",
    width = width,

    selectInput(ns("dge_method"), label = "Method", choices = NULL),
    pickerInput(
      ns("features"),
      label = "Feature Class",
      choices = NULL,
      multiple = TRUE,
      options = list(
        `selected-text-format`= "count",
        `count-selected-text` = "{0} classes chosen"
      )),
    numericInput(ns("treatfc"), label = "Min. Fold Change",
                 value = 1, min = 1, max = 10, step = 0.25),
    prettyCheckbox(ns("sample_weights"), label = "Use Sample Weights",
                   status = "primary")
    # , tags$h4("Filter Strategy"),
    # wellPanel(
    #   fdgeFeatureFilterUI(ns("ffilter"), ..., debug = debug))
  )
}

#' @noRd
initialized.FdgeRunOptions <- function(x, ...) {
  # Waiting on: FacileAnalysis/issues/8
  dge.method <- x$dge_method()
  !unselected(dge.method) &&
    dge.method != "ranks" &&
    nrow(x$feature_filter$universe()) > 0
}

# Filtering Options ============================================================

#' @noRd
#' @param assay_mod [FacileShine::assaySelect()] module, ie. an
#'   `AssaySelectInput` object.
fdgeFeatureFilter <- function(input, output, session, rfds, model,
                              assay_module, ..., debug = FALSE) {

  # TODO: Support returning a vector of feature id's
  filter_method <- reactive("default")

  min_count <- reactive({
    input$min_count
  })
  min_total_count <- reactive({
    input$min_total_count
  })
  min_expr <- reactive({
    input$min_expr
  })

  # Show/hide the appropropriate options. I would have used
  # shinyjs::toggleElement, but that's not working inside the dropdown
  # this is embedded in
  # observe({
  #   toggleElement("countopts", condition = count_options())
  #   toggleElement("expropts", condition = !count_options())
  # })
  output$show_count_options <- reactive({
    assay <- assay_module$assay_info()$assay_type
    if (assay %in% c("rnaseq", "umi", "isoseq")) "yes" else "no"
  })
  outputOptions(output, "show_count_options", suspendWhenHidden = FALSE)

  vals <- list(
    method = filter_method,
    min_count = min_count,
    min_total_count = min_total_count,
    min_expr = min_expr)
}

#' @noRd
#' @importFrom shiny conditionalPanel NS numericInput tagList tags
fdgeFeatureFilterUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    conditionalPanel(
      condition = "output.show_count_options == 'yes'", ns = ns,
      tags$div(id = ns("countopts"),
               numericInput(ns("min_count"), label = "Min. Count",
                            value = 10, min = 1, max = 50, step = 5),
               numericInput(ns("min_total_count"), label = "Min. Total Count",
                            value = 15, min = 1, max = 100, step = 10))),
    conditionalPanel(
      condition = "output.show_count_options == 'no'", ns = ns,
      tags$div(id = ns("expropts"),
               numericInput(ns("min_expr"), label = "Min. Expression",
                            value = 1, min = 0.25, max = 100, step = 5))))
}
