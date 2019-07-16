#' This module will retun a list of arguments to configure the ffsea method
#' for the specific subclass of an AnalysisResult that `ares` belongs to.
#' The list that is returned here will not include `ares` itself.
#'
#' @noRd
#' @export
#' @importFrom shiny renderUI
#' @param ares a reactive that contains a FacileAnalysisResult
#' @return An `FfseaRunOptions` list, where `$args()` is a reactive list, with
#'   name=value pairs set to the arguments to run the appropriate `ffsea.*`
#'   method for `aresult`.
ffseaRunOpts <- function(input, output, session, rfds, aresult, ...,
                         debug = FALSE) {
  ns <- session$ns
  state <- reactiveValues(
    aclass = "__initializing__",
    default_args = "__initializing__")

  ares <- reactive({
    req(initialized(aresult))
    faro(aresult)
  })

  observeEvent(ares(), {
    ares. <- req(ares())
    aclass <- class(ares.)[1L]

    # always reset to default params when a new analysis result is provided,
    # but we may want to change this later.
    if (TRUE || state$aclass != aclass) {
      rm.ui <- names(state$default_args)
      state$aclass <- aclass
      uistuff <- renderFfseaRunOptsUI(ares., ns)
      state$default_args <- uistuff$args
      # renderUI("ui", uistuff$ui)
      output$ui <- renderUI(uistuff$ui)
    }
  }, priority = 5)

  # extracts the values from the UI for this object
  args <- reactive({
    dargs <- state$default_args
    req(!unselected(dargs))
    lapply(names(dargs), function(name) input[[name]])
  })

  # 1. remove the UI already there
  # 2. create the new one
  vals <- list(
    args = args,
    aclass = reactive(state$aclass),
    .state = state,
    .ns = session$ns)
  class(vals) <- "FfseaRunOptions"
  vals
}

#' Creates the UI panel for the FacileAnalysisResult-specific ffsea options.
#'
#' This function creates the UI `tagList` for the options required to run
#' `ffsea` over a specific analysis result. The values are initialized to
#' their defaults. The `tagList` and list of `args` are returned.
#'
#' @noRd
#' @export
#' @importFrom shiny NS tags uiOutput
#' @importFrom shinyWidgets dropdownButton
#' @return a list with `$ui` for the tagList of interface components and
#'   `$args`, which is a list of name/value pairs for the default arguments
#'   of the ffsea.* function implementation.
ffseaRunOptsUI <- function(id, width = "300px", ..., debug = FALSE) {
  ns <- NS(id)
  dropdownButton(
    inputId = ns("opts"),
    icon = icon("sliders"),
    status = "primary", circle = FALSE,
    width = width,
    tags$div(
      id = ns("ffseaRunOptsContainer"),
      # style = "height: 400px",
      uiOutput(ns("ui"))))
}

# Helper Functions =============================================================

#' @noRd
#' @export
initialized.FfseaRunOptions <- function(x, ...) {
  args <- x$args
  is.list(args) && length(args) > 0L && !any(sapply(args, is.na))
}

#' @noRd
#' @export
#' @importFrom shiny tagList
#' @param x the FacileAnalysisResult to buil a UI for
#' @param ns the namespace object to use for ui id generation
renderFfseaRunOptsUI <- function(x, ns, ...) {
  UseMethod("renderFfseaRunOptsUI", x)
}

#' @noRd
#' @export
#' @importFrom shiny checkboxInput numericInput selectInput tagList tags
renderFfseaRunOptsUI.FacileTtestAnalysisResult <- function(x, ns, ...) {
  args <- default_ffsea_args(x)
  # We can rank by either logFC or t-statistics (if they are available)
  rank.opts <- intersect(c("logFC", "t"), colnames(tidy(x)))

  # NOTE: We may want to update this and even have the UI update with the number
  # of features that pass these filtering criteria, and what have you.

  ui <- tagList(
    tags$h4("Enrichment Options"),
    numericInput(ns("min_logFC"), "Min logFC", value = args$min_logFC,
                 min = 0, max = 10, step = 0.5),
    numericInput(ns("max_padj"), "Max FDR", value = args$max_padj,
                 min = 0, max = 1, step = 0.05),
    tags$h4("Preranked Options"),
    selectInput(ns("rank_by"), "Rank By", choices = rank.opts,
                selected = rank.opts[1]),
    checkboxInput(ns("signed"), "Signed Ranks", value = args$signed))

  list(args = args, ui = ui)
}

#' @noRd
#' @export
#' @importFrom shiny tagList tags
renderFfseaRunOptsUI.FacileAnovaAnalysisResult <- function(x, ns, ...) {
  args <- default_ffsea_args(x)
  ui <- tagList()
  list(args = args, ui = ui)
}

#' @noRd
#' @export
#' @importFrom shiny tagList tags
renderFfseaRunOptsUI.FacilePcaAnalysisResult <- function(x, ns, ...) {
  args <- default_ffsea_args(x)
  ui <- tagList()
  list(args = args, ui = ui)
}

