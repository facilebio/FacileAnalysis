#' Wraps FacileAnalysis modules into a gadget.
#'
#' If `x` is a `facile_frame`, than the analysis is strictly performed over the
#' samples defined by the `facile_frame`. If it is a `FacileDataStore`, then
#' a UI is provided to filter into the sample space by the
#' [FacileShine::filteredReactiveFacileDataStoreUI()] module.
#'
#' @importFrom miniUI
#'   gadgetTitleBar
#'   miniPage
#' @importFrom multiGSEA failWith
#' @importFrom FacileShine
#'   filteredReactiveFacileDataStore
#'   filteredReactiveFacileDataStoreUI
#'   ReactiveFacileDataStore
#' @importFrom shiny
#'   browserViewer
#'   callModule
#'   dialogViewer
#'   runGadget
#'   tagList
#'   tags
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyWidgets useSweetAlert
#' @param analysisModule the server function for the gadget
#' @param analysisUI the UI function for the gadget
#' @param x the FacileDataStore or a facile sample descriptor from one.
#' @param user the username/id for the analysis
#' @param title the title of the gadget
#' @param height,width height and width of gadget, defaults to 600 and 1000,
#'   respectively.
#' @param viewer `"dialog"` (default), `"browser"`, or `"pane"`. See description
#'   from [shiny::runGadget]
#' @return a subclass of `FacileAnalysisGadgetResult` which holds the result
#'   from the `analysisModule`, as well as am `annotation` tibble, that can hold
#'   information from brushing of smaples (or features) -- depending on what
#'   the gadget is working over.
frunGadget <- function(analysisModule, analysisUI, x, user = Sys.getenv("USER"),
                       title = "Facile Analysis Gadget",
                       height = 600, width = 1000, viewer = "pane", ...,
                       debug = FALSE) {
  bs4dash <- getOption("facile.bs4dash")
  options(facile.bs4dash = FALSE)
  on.exit(options(facile.bs4dash = bs4dash))

  assert_function(analysisModule)
  assert_function(analysisUI)

  viewer <- gadget_viewer(viewer, title, width, height)

  if (is(x, "facile_frame")) {
    fds. <- fds(x)
    samples. <- x
    sample.filter <- FALSE
  } else {
    sample.filter <- TRUE
    fds. <- x
    samples. <- samples(x) %>% collect(n = Inf)
  }

  assert_class(fds., "FacileDataStore")
  assert_class(samples., "facile_frame")

  if (sample.filter) {
    ui.content <- tagList(
      filteredReactiveFacileDataStoreUI("ds"),
      tags$hr(),
      analysisUI("analysis", ..., debug = debug))
  } else {
    ui.content <- analysisUI("analysis", ..., debug = debug)
  }

  ui <- miniPage(
    useShinyjs(),
    useSweetAlert(),
    gadgetTitleBar(title),
    miniContentPanel(ui.content),
    NULL)

  server <- function(input, output, session) {
    rfds <- ReactiveFacileDataStore(fds., "ds", user = user, samples = samples.)
    analysis <- callModule(analysisModule, "analysis", rfds, ..., debug = debug)

    observeEvent(input$done, {
      annotation <- FacileShine:::.empty_feature_annotation_tbl()
      result. <- failWith(list(), unreact(result(analysis)))
      result.[["gadget"]] <- list(
        annotation = annotation)
      class(result.) <- classify_as_gadget(result.)
      stopApp(invisible(result.))
    })
    observeEvent(input$cancel, {
      stopApp(invisible(NULL))
    })
  }

  runGadget(ui, server, viewer = viewer, stopOnCancel = FALSE)
}


#' @noRd
#' @importFrom shiny browserViewer dialogViewer paneViewer
gadget_viewer <- function(type = c("dialog", "browser", "pane"),
                          title = "Facile Gadget",
                          width = 800, height = 600, minHeight = height,
                          ...) {
  type <- match.arg(type)
  switch(type,
         dialog = dialogViewer(title, width, height),
         browser = browserViewer(),
         pane = paneViewer(minHeight = minHeight))
}
