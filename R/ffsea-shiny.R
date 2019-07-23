#' Full interactive executation of GSEA, soup to nuts.
#'
#' For now, feature set enrichment analysis are only performed downstream of
#' a `FacileAnalysisResult`. This requires that we have an ffsea method defined
#' for the specific class of the result, ie. `ffsea.FacilePcaAnalysisResult`.
#'
#' @rdname interactive-ffsea
#'
#' @export
#' @param x A FacileAnalysisResult that has an implemented `ffsea.*` method
#' @param gdb A `GeneSetDb` to use for the FSEA.
#' @examples
#' gdb <- multiGSEA::exampleGeneSetDb()
#' dge.crc <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor", denom = "normal", fixed = "sex") %>%
#'   fdge(method = "voom")
#' if (interactive()) {
#'   ffseaGadget(dge.crc, gdb)
#' }
ffseaGadget <- function(x, gdb, title = "Feature Set Enrichment Analysis",
                        height = 800, width = 1000, ..., debug = FALSE) {
  assert_class(x, "FacileAnalysisResult")
  assert_class(gdb, "GeneSetDb")
  frunGadget(ffseaAnalysis, ffseaAnalysisUI, x, aresult = x, gdb = gdb,
             title = title, height = height, width = width, ...,
             retval = "faro", debug = debug)
}

#' A moodule that encapsulates configuring and running ffsea, and a view to
#' interact with the results.
#'
#' @noRd
#' @export
ffseaAnalysis <- function(input, output, session, rfds, aresult, gdb = NULL,
                          ..., debug = FALSE) {
  res <- callModule(ffseaRun, "run", rfds, aresult, gdb, ..., debug = debug)
  view <- callModule(ffseaView, "view", rfds, res, ..., debug = FALSE)

  # Only show the view UI when there is an FfseaAnalysisResult ready
  observe({
    toggleElement("viewbox", condition = initialized(res))
  })

  vals <- list(
    main = res,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacileFseaAnalysisResultContainer",
                   "ReactiveFacileAnalysisResultContainer")
  vals
}

#' @noRd
#' @export
#' @importFrom shinyjs hidden
ffseaAnalysisUI <- function(id, ...) {
  ns <- NS(id)

  tagList(
    tags$div(
      id = ns("runbox"),
      box(title = "Configure Feature Set Analysis", width = 12,
          ffseaRunUI(ns("run")))),
    hidden(
      tags$div(
        id = ns("viewbox"),
        box(title = "Feature Set Analysis Results", solidHeader = TRUE,
            width = 12, ffseaViewUI(ns("view"), debug = debug)))))
}

# Run ==========================================================================

#' @section Custom Run Configuration:
#' The options presented to the user for running a feature set enrichment
#' analysis on a `FacilePcaAnalysisResult` will be different than the ones
#' made available for an enrichment analysis over a `FacileTtestAnalysisResult`
#' or even a `FacileAnovaAnalysisResult`.
#'
#' As such, each type of result should define a UI that accepts the appropriate
#' parameters for its corresponding `ffsea.*` method, and a server function
#' that extract and invokes the function.
#'
#' @rdname interactive-ffsea
#' @export
#' @importFrom shiny eventReactive withProgress
#' @importFrom shinyWidgets updatePickerInput
#'
#' @param aresult A `FacileAnalysisResult` that has a `ffsea.*` method defined.
#' @param gdb A GeneSetDb object to use for testing.
ffseaRun <- function(input, output, session, rfds, aresult, gdb = NULL, ...,
                     debug = FALSE) {
  ares <- reactive({
    req(initialized(aresult))
    faro(aresult)
  })

  # Enable user to configure the GeneSetDb used for testing
  gdb. <- callModule(geneSetDbConfig, "gdb", rfds, aresult = aresult,
                     gdb = gdb, ..., debug = debug)

  # When the AnalysisResult changes, update the runopts UI based on the specific
  # subclass of the AnalysisResult we have at play
  runopts <- callModule(ffseaRunOpts, "runopts", rfds, aresult = aresult, ...,
                        debug = debug)


  # Updates the set-enrichment methods when the analysis result changes.
  available_methods <- reactive({
    ares. <- req(ares())
    ffsea_methods(ares.)
  })

  observeEvent(available_methods(), {
    methods <- req(available_methods())
    choices <- split(methods[["method"]], methods[["type"]])
    selected <- unname(sapply(choices, "[", 1L))
    opts <- NULL
    updatePickerInput(session, "ffsea_methods", selected = selected,
                      choices = choices, choicesOpt = opts)
  })

  runnable <- reactive({
    !unselected(input$ffsea_methods) &&
      initialized(ares()) &&
      initialized(gdb.) &&
      initialized(runopts)
  })

  observe({
    runnable. <- runnable()
    ftrace("runnable: ", as.character(runnable.))
    toggleState("runbtn", condition = runnable.)
  })

  observe({
    ftrace("Run button pressed: ", as.character(input$runbtn))
  })

  fsea_res <- eventReactive(input$runbtn, {
    req(runnable())
    args. <- runopts$args()
    gdb. <- gdb.$gdb()
    args <- c(list(x = ares(), gdb = gdb.), args.)
    withProgress({
      do.call(ffsea, args)
    }, message = "Running Enrichment Analysis")
  })

  vals <- list(
    faro = fsea_res,
    .ns = session$ns)

  # TODO: fix class hierarchy
  classes <- c("ReactiveFacileFseaAnalysisResult",
               "ReactiveFacileAnalysisResult",
               "FacileFseaAnalysisResult")

  class(vals) <- classes

  vals
}

#' @noRd
#' @export
#' @importFrom shinyWidgets pickerInput
#' @importFrom shiny actionButton column fluidRow tags tagList
ffseaRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  tagList(
    tags$div(
      id = ns("gdb-container"),
      geneSetDbConfigUI(ns("gdb"))),
    fluidRow(
      column(
        3,
        pickerInput(ns("ffsea_methods"), "Analysis Methods", choices = NULL,
                    multiple = TRUE)),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          ffseaRunOptsUI(ns("runopts"), width = "300px"))),
      column(1, actionButton(ns("runbtn"), "Run", style = "margin-top: 1.7em"))
    )
  )
}

# View =========================================================================

#' Responsible for the shiny view of a FacileFseaAnalysisResult
#'
#' @noRd
#' @export
#' @importFrom multiGSEA failWith
#' @importFrom multiGSEA.shiny
#'   geneSetContrastView
#'   mgGeneSetSummaryByGene
#'   mgResultFilter
#'   mgTableBrowser
#'   MultiGSEAResultContainer
#'   summaryHTMLTable.multiGSEA
#'   updateActiveGeneSetInContrastView
#' @importFrom shiny validate
#' @param rfds the reactive facile data store
#' @param ares The `FacileFseaAnalysisResult`
ffseaView <- function(input, output, session, rfds, aresult, ...,
                      debug = FALSE) {
  state <- reactiveValues(
    gsview_select = tibble(assay_name = character(), feature_id = character()),
    set_select = tibble(collection = character(), name = character())
  )

  ares <- reactive({
    req(initialized(aresult))
    faro(aresult)
  })

  mgc <- reactive({
    res <- req(ares())
    mgres <- result(res)
    # TODO: validate() isn't working here, only works inside a
    # `output$xxx <- render*({})` block, which is generating outpout into the
    # shiny app
    validate(
      need(is(mgres, "MultiGSEAResult"), "MultiGSEAResult can't be found")
    )

    MultiGSEAResultContainer(mgres)
  })

  gs_result_filter <- callModule(mgResultFilter, "mg_result_filter", mgc)

  # Overview Tab ...............................................................
  output$gseaMethodSummary <- renderUI({
    mgc. <- req(mgc())
    tagList(
      tags$h4("GSEA Analyses Overview"),
      summaryHTMLTable.multiGSEA(mgc.$mg, mgc.$methods,
                                 gs_result_filter()$fdr(),
                                 p.col = "padj.by.collection")
    )
  })

  # GSEA Results Tab ...........................................................
  gs_viewer <- callModule(geneSetContrastView, "geneset_viewer",
                          mgc, maxOptions=500, server=TRUE)

  # A table of GSEA statistics/results for the given method and fdr threshold
  # The table is wired to the gs_viewer so that row clicks can signal updates
  # to the contrast viewer
  gs_table_browser <- callModule(mgTableBrowser, "mg_table_browser", mgc,
                                 method=gs_result_filter()$method,
                                 fdr=gs_result_filter()$fdr,
                                 server=TRUE)
  # clicks on gsea result table update the contrast view
  observeEvent(gs_table_browser$selected(), {
    .mgc <- req(mgc())
    geneset <- req(gs_table_browser$selected())
    updateActiveGeneSetInContrastView(session, gs_viewer, geneset, .mgc)
  })

  # A table of other genesets that brushed genes in the contrast viewer
  # belong to. This table is also wired to the contrast viewer, so that
  # a click on a row of the table will update the contrast view, too.
  other_genesets_gsea <- callModule(mgGeneSetSummaryByGene,
                                    "other_genesets_gsea",
                                    mgc, features = gs_viewer()$selected,
                                    method = gs_result_filter()$method,
                                    fdr = gs_result_filter()$fdr)

  vals <- list(
    selected_features = reactive(state$gsview_select),
    selected_sets = reactive(state$set_select),
    .ns = session$ns)

  vals
}

#' @noRd
#' @export
#' @importFrom shiny fluidRow NS tags uiOutput wellPanel
#' @importFrom multiGSEA.shiny
#'   geneSetContrastViewUI
#'   mgGeneSetSummaryByGeneUI
#'   mgResultFilterUI
#'   mgTableBrowserUI
ffseaViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  tagList(
    wellPanel(mgResultFilterUI(ns("mg_result_filter"))),

    tags$div(
      style="margin-bottom: 10px; padding: 5px; background-color: white",
      title='GSEA Results',
      uiOutput(ns("gseaMethodSummary"))),

    fluidRow(
      column(
        5, style="padding: 0",
        wellPanel(geneSetContrastViewUI(ns("geneset_viewer")))),
      column(
        7, mgTableBrowserUI(ns("mg_table_browser")))),

    fluidRow(
      column(
        12,
        tags$h4("Other Gene Sets with Selected Genes"),
        mgGeneSetSummaryByGeneUI(ns("other_genesets_gsea"))))
  )
}
