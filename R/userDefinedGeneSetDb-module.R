#' Creates a GeneSetDb from a user-specfied gene set definition
#'
#' This module provides an upload button that allows a user to upload a
#' table of gene set definitions in [mutiGSEA::GeneSetDb()] format. Minimal
#' validation checks are implemented.
#'
#' @export
#' @importFrom multiGSEA GeneSetDb
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
userDefinedGeneSetDb <- function(input, output, session, ...) {

  empty.def <- tibble(
    collection = character(), name = character(), featureId = character())

  state <- reactiveValues(
    dat = empty.def)

  observeEvent(input$upload, {
    type <- input$upload$type
    path <- input$upload$datapath
    ext <- file_ext(path)

    validate(
      need(ext %in% c("csv", "xlsx"), "Only csv or xlsx files are supported")
    )

    if (ext == "csv") {
      dat <- tryCatch(
        read.csv(path, stringsAsFactors = FALSE),
        error = function(e) NULL)
    } else {
      dat <- tryCatch(
        read_excel(path),
        error = function(e) NULL)
    }

    validate(
      need(is.data.frame(dat), "Error parsing geneset definition file")
    )

    req.cols <- c("collection", "name", "featureId")
    missed <- setdiff(req.cols, colnames(dat))
    validate(
      need(
        length(missed) == 0L,
        sprintf("Missing columns: %s", paste(missed, collapse = ","))))

    dat[["featureId"]] <- as.character(dat[["featureId"]])
    state$dat <- dat %>%
      mutate(collection = as.character(collection),
             name = as.character(name),
             featureId = as.character(featureId))
  })

  gdb <- reactive({
    dat. <- state$dat
    out <- if (is.null(dat.) || nrow(dat.) == 0) NULL else GeneSetDb(dat.)
    out
  })

  observe({
    toggleState("rm", condition = is(gdb(), "GeneSetDb"))
  })

  observeEvent(input$rm, {
    state$dat <- empty.def
  })

  vals <- list(
    gdb = gdb,
    .ns = session$ns)
  class(vals) <- "ReactiveGeneSetDb"
  vals
}

#' @noRd
#' @export
#' @importFrom shiny fileInput tagList
userDefinedGeneSetDbUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  inline.style <- "display: inline-block; vertical-align:top;"

  tagList(
    tags$div(
      style = inline.style,
      fileInput(
        ns("upload"), "Additional Gene Set Definition File",
        multiple = FALSE,
        accept = c(
          # CSV stuff
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv",
          # xlsx stuff (not .xls)
          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          ".xlsx"))),
    tags$div(
      style = paste(inline.style, "margin-top: 1.8em"),
      actionButton(ns("rm"), "Reset")))

}
