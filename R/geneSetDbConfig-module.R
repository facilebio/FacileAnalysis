#' Module to instantiate and configure a GeneSetDb for a FacileAnalysisResult
#'
#' Currently we accept
#' Need to provide user with the collection,geneset combinations that are
#' appropriate for the features in the assay_tyupe of `ares`.
#'
#' TODO: There's a lot to do in order to make a shiny interaface to a GeneSetDb
#'
#' @noRd
#' @export
#' @importFrom multiGSEA append conform GeneSetDb featureIds geneSets
#' @importFrom shiny callModule need reactive renderUI validate
#' @param rfds `ReactiveFacileDataSet`
#' @param aresult A `FacileAnalysisResult` to run through ffsesa
#' @param gdb An optional "super" GeneSetDb that we will subset from.
geneSetDbConfig <- function(input, output, session, rfds, aresult, gdb = NULL,
                            ..., debug = FALSE) {

  ares <- reactive({
    req(initialized(aresult))
    faro(aresult)
  })

  # features used in the analysis
  afeatures <- reactive({
    req(ares()) %>%
      features()
  })

  # The default "superset" GeneSetDb. This is passed into modules to run
  # ffsea. In the future we need to plug into a GeneSetDb provider to
  # define this.
  gdb.base <- reactive({
    if (is(gdb, "GeneSetDb")) {
      gdb. <- gdb
    } else if (is(gdb, "ReactiveGeneSetDb")) {
      gdb. <- gdb$gdb()
    } else if (is.null(gdb)) {
      gdb. <- NULL
    } else {
      stop("We need to connect to some remote GeneSetDb provider ...")
    }
    gdb.
  })

  gdb.user <- callModule(userDefinedGeneSetDb, "user_gdb")

  gdb.go <- reactive({
    finfo <- req(afeatures())
    base. <- gdb.base()
    user. <- gdb.user$gdb()

    if (is.null(base.)) {
      validate(
        need(is(user., "GeneSetDb"),
             "No Default GeneSetDb provided, please upload custom gene sets"))
      out <- user.
    } else {
      out <- base.
      if (!is.null(user.)) {
        out <- append(out, user.)
      }
    }

    if (!is.null(out)) {
      out <- conform(out, finfo[["feature_id"]])
    }
    out
  })

  output$message <- renderUI({
    gdb. <- gdb.go()
    if (!is(gdb., "GeneSetDb")) {
      msg <- "No GeneSetDb defined"
    } else {
      gs <- geneSets(gdb.)
      msg <- glue(
        "GeneSetDb defined: {ngs} gene sets, {ncol} collections, {ngene} genes",
        ngs = nrow(gs),
        ncol = length(unique(gs[["collection"]])),
        ngene = length(featureIds(gdb.)))
    }
    tags$p(msg)
  })

  vals <- list(
    gdb = gdb.go,
    .ns = session$ns)
  class(vals) <- c("ReactiveGeneSetDb")

  vals
}

#' @noRd
#' @export
#' @importFrom shiny NS tags uiOutput
geneSetDbConfigUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    tags$div(id = ns("userdb"), userDefinedGeneSetDbUI(ns("user_gdb"))),
    tags$div(id = ns("msg"), uiOutput(ns("message"))))
}

#' @noRd
#' @export
initialized.ReactiveGeneSetDb <- function(x, ...) {
  gdb. <- x$gdb()
  # is(gdb., "GeneSetDb") && is.conformed(gdb.)
  is(gdb., "GeneSetDb")
}

