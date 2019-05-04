#' Takes care of the tedium of wrapping a facilemodule to run as a gadget.
#'
#' We will implement barebones analysis gadgets that delegate down into here.
#'
#' @importFrom miniUI
#'   gadgetTitleBar
#'   miniPage
#' @importFrom multiGSEA failWith
#' @importFrom FacileShine ReactiveFacileDataStore
#' @importFrom shiny
#'   browserViewer
#'   callModule
#'   dialogViewer
#'   runGadget
frunGadget <- function(analysisModule, analysisUI, x, user = Sys.getenv("USER"),
                       title = "Facile Analysis Gadget",
                       height = 600, width = 1000, viewer = "dialog", ...,
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
  } else {
    fds. <- x
    samples. <- samples(x)
  }

  assert_class(fds., "FacileDataStore")
  assert_class(samples., "facile_frame")

  ui <- miniPage(
    gadgetTitleBar(title),
    miniContentPanel(analysisUI("analysis", ..., debug = debug)),
    NULL)

  server <- function(input, output, session) {
    rfds <- ReactiveFacileDataStore(fds., "ds", user = user, samples = samples.)
    analysis <- callModule(analysisModule, "analysis", rfds, ..., debug = debug)

    observeEvent(input$done, {
      annotation <- FacileShine:::.empty_feature_annotation_tbl()
      result. <- failWith(list(), unreact(analysis$result$result()))
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
                          width = 800, height = 600, minHeight = height) {
  type <- match.arg(type)
  switch(type,
         dialog = dialogViewer(title, width, height),
         browser = browserViewer(),
         pane = paneViewer(minHeight = minHeight))
}
