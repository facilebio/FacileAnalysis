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
#' @examples
#'
#' gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
#' dge.ttest <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor", denom = "normal", fixed = "sex") %>%
#'   fdge(method = "voom")
#' gsea.ttest <- ffsea(dge.ttest, gdb, methods = c("cameraPR", "fgsea"))
#' if (interactive()) {
#'   shine(ttest.gsea)
#' }
ffseaGadget <- function(x, title = "Feature Set Enrichment Analysis",
                        height = 800, width = 1000, ...) {
  assert_class(x, "FacileAnalysisResult")
  frunGadget(ffseaAnalysis, ffseaAnalysisUI, x, title = title,
             height = height, width = width, ...)
}

#' A moodule that encapsulates configuring and running ffsea, and a view to
#' interact with the results.
#'
#' @noRd
#' @export
ffseaAnalysis <- function(input, output, session, rfds, aresult, ...,
                          debug = FALSE) {
  res <- callModule(ffseaRun, "run", rfds, aresult, ..., debug = debug)
  view <- callModule(ffseaView, "view", rfds, res, ..., debug = FALSE)

  # Only show the view UI when there is an FfseaAnalysisResult ready
  observe({
    res. <- req(faro(res))
    show <- is(res., "FacileFseaAnalysisResult")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    main = res,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacileFseaAnalysisResultContainer",
                   "ReactiveFacileAnalysisResultContainer")
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
          ffseaRunUI("run"))),
    hidden(
      tags$div(
        id = ns("viewbox"),
        box(title = NULL, solidHeader = TRUE, width = 12,
            ffseaViewUI(ns("view"), debug = debug)))))
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
#'
#' @param aresult A `FacileAnalysisResult` that has a `ffsea.*` method defined.
#' @param gdb A GeneSetDb object to use for testing.
#' @export
ffseaRun <- function(input, output, session, rfds, aresult, gdb = NULL, ...,
                     debug = FALSE) {

  ares <- reactive({
    req(initialized(ares))
    faro(ares)
  })

  # When the AnalysisResult changes, update the runopts UI based on the specific
  # subclass of the AnalysisResult we have at play
  runopts <- callModule(ffseaRunOpts, "runopts", rfds, ares, ...)

  gdb. <- callModule(geneSetDbConfig, "gdb", rfds, ares, ...)

  # When the FacileAnalysisResult changes, we need to update the types of
  # methods that are available for analysis on this thing.
  observe({
    ares. <- req(ares)

  })

  fsea_methods <- reactive({
    methods <- input$ffsea_methods
  })

  runnable <- reactive({
    initialized(ares()) &&
      initialized(runopts) &&
      length(fsea_methods()) > 0
  })

  observe({
    toggleState("run", condition = runnable())
  })

  fseares <- eventReactive(input$run, {
    req(runnable())
    args. <- runopts$args()
    gdb. <- gdb$gdb()
    args <- c(list(x = ares(), gdb = gdb$gdb()), args.)
    out <- do.call(ffsea, args)
  })


  vals <- list(
    faro = fseares,
    .ns = session$ns)

  # TODO: fix class hierarchy
  classes <- c("ReactiveFacileFseaAnalysisResult",
               "FacileFseaAnalysisResult",
               "ReactiveFacileAnalysisResult")

  class(vals) <- classes
  vals
}

#' @noRd
#' @export
#' @importFrom shinyWidgets pickerInput
ffseaRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  tagList(
    tags$div(id = "gdb-container", geneSetDbConfigUI(ns("gdb"))),
    fluidRow(
      column(
        3,
        pickerInput(ns("ffsea_methods"), "Analysis Methods", choices = NULL)),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          ffseaRunOptsUI(ns("runopts"), width = "300px"))),
      column(1, actionButton(ns("run"), "Run", style = "margin-top: 1.7em"))
    )
  )
}

#' @noRd
#' @export
initialized.FfseaRunOptions <- function(x, ...) {
  test_list(x, names = "unique") && all(sapply(vals, is, "reactive"))
}

#' @noRd
#' @export
ffseaRunOpts <- function(x, ...) {
  UseMethod("ffseaRunOpts", x)
}

ffseaRunOpts.default <- function(x, ...) {
  warning("ffseaRunOpts not implemented for: ", class(x)[1L])
  list(arg1 = NULL, arg2 = NULL)
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
ffseaView <- function(input, output, session, rfds, ares, ...,
                      debug = FALSE) {
  state <- reactiveValues(
    gsview_select = tibble(assay_name = character(), feature_id = character()),
    set_select = tibble(collection = character(), name = character())
  )

  # assert_class()
  fsea_res <- reactive({
    req(initialized(ares))
    faro(ares)
  })

  mgc <- reactive({
    res <- req(fsea_res())
    mgres <- result(res)

    validate(
      need(is(mgres, "MultiGSEAResult"), "MultiGSEAResult can't be found")
    )

    MultiGSEAResultContainer(mgres)
  })

  gs_result_filter <- callModule(mgResultFilter, "mg_result_filter", mgc)

  # Overview Tab ===============================================================
  output$gseaMethodSummary <- renderUI({
    mgc. <- req(mgc())
    tagList(
      tags$h4("GSEA Analyses Overview"),
      summaryHTMLTable.multiGSEA(mgc.$mg, mgc.$methods,
                                 gs_result_filter()$fdr(),
                                 p.col = "padj.by.collection")
    )
  })

  # GSEA Results Tab ===========================================================
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
