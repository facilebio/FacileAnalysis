#' Full interactive executation of GSEA, soup to nuts.
#'
#' I think this has to be loaded with an already run FacileAnalysisResult
#' that we have an ffsea method for, ie. configure a GSEA from an already baked
#' DGE or PCA AnalysisResult.
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

}

#' @noRd
#' @export
ffseaAnalysisUI <- function(id, ...) {
  ns <- NS(id)
}

#' @noRd
#' @export
ffseaRun <- function(input, output, session, rfds, ares, ..., debug = FALSE) {

  vals <- list(
    faro = reactive(NULL),
    .ns = session$ns)

  classes <- c("FacileFseaAnalysisResult", "...more...")
  class(vals) <- classes
  vals
}

#' @noRd
#' @export
ffseaRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
}

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
