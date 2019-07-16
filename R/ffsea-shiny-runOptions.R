#' This module will retun a list of arguments to configure the ffsea method
#' for the specific subclass of an AnalysisResult that `ares` belongs to.
#' The list that is returned here will not include `ares` itself.
#'
#' @noRd
#' @export
#' @importFrom shiny renderUI
#' @param ares a reactive that contains a FacileAnalysisResult
ffseaRunOpts <- function(input, output, session, rfds, ares, ...,
                         debug = FALSE) {
  ns <- session$ns
  state <- reactiveValues(
    aclass = "__initializing__",
    default_args = "__initializing__",
  )

  observeEvent(ares(), {
    ares. <- req(ares())
    aclass <- class(ares.)[1L]

    if (state$aclass != aclass) {
      rm.ui <- names(state$default_args)
      state$aclass <- aclass
      uistuff <- renderFfseaRunOptsUI(ares., ns)
      state$default_args <- uistuff$args
      renderUI("ui", uistuff$ui)
    }
  }, priority = 5)

  observeEvent(state$aclass, {
    # Add new UI elements

  })

  # ares is already a reactive FacileAnalysisResult
  # When it changes, we have to update the GUI and get its values

  argnames <- reactive(names(state$default_args))

  # extracts the values from the UI for this object
  args <- reactive({
    lapply(argnames(), function(name) input[[name]])
  })

  # 1. remove the UI already there
  # 2. create the new one
  vals <- list(
    args = args,
    aclass = reactive(state$aclass),
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
      uiOutput(ns("ui"))))
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
#' @importFrom shiny tagList tags
renderFfseaRunOptsUI.FacileTtestAnalysisResult <- function(x, ns, ...) {
  args <- default_ffsea_args(x)
  ui <- tagList()
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

