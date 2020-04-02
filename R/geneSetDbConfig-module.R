#' Module to instantiate and configure a GeneSetDb for a FacileAnalysisResult
#'
#' This accepts a "supserset-of-a-GeneSetDb" that can be added to by a user
#' defined set of gene sets via the [multiGSEA.shiny::userDefinedGeneSetDb()]
#' module. If this module is not initialized with such a `GeneSetDb` via the
#' `gdb` argument, then the only `GeneSetDb` available will be the one that is
#' uploaded by the user.
#'
#' Need to provide user with the collection,geneset combinations that are
#' appropriate for the features in the assay_tyupe of `aresult`. Although we are
#' currently useing `aresult` to decipher the feature space of the GeneSetDb
#' to generate, this should probably as it is too specific to Facile, and this
#' module more appropriately belongs in the multiGSEA.shiny package.
#'
#' @section Refactoring to support connection to a GeneSetDb provider:
#' We should be able to query a GeneSetDb provider to ask for a catlog of
#' gene set collections by providing the organism and feature type that we want
#' gene sets for.
#'
#' @noRd
#' @export
#' @importFrom multiGSEA append conform GeneSetDb featureIds geneSets
#' @importFrom multiGSEA.shiny userDefinedGeneSetDb
#' @importFrom shiny callModule need reactive renderUI validate
#'
#' @param rfds `ReactiveFacileDataSet`
#' @param aresult A `FacileAnalysisResult` to run through ffsesa
#' @param gdb An optional "super" GeneSetDb that we will subset from.
#' @return A ReactiveGeneSetDb (list) with a configured GeneSetDb ready for
#'   analysis over `aresult` in `$gdb()`
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
    req(ares())
    if (is(gdb, "reactive")) {
      gdb. <- gdb()
      if (!is(gdb., "GeneSetDb")) {
        gdb. <- NULL
      }
    } else if (is(gdb, "GeneSetDb")) {
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

  gdb.user <- callModule(userDefinedGeneSetDb, "user_gdb", ...)

  gdb.go <- reactive({
    finfo <- req(afeatures())
    base.gdb <- gdb.base()
    user.gdb <- gdb.user$gdb()

    if (is.null(base.gdb)) {
      validate(
        need(is(user.gdb, "GeneSetDb"),
             "No Default GeneSetDb provided, please upload custom gene sets"))
      out <- user.gdb
    } else {
      out <- base.gdb
      if (!is.null(user.gdb)) {
        out <- append(out, user.gdb)
      }
    }

    # if (!is.null(out)) {
    #   out <- conform(out, finfo[["feature_id"]])
    # }
    out
  })

  # TODO: Implemenet a better messaging system for the geneSetDbConfig here.
  output$message <- renderUI({
    gdb. <- gdb.go()
    if (!is(gdb., "GeneSetDb")) {
      msg <- "No GeneSetDb defined"
    } else {
      gs.all <- geneSets(gdb., active.only = FALSE)
      gs.active <- filter(gs.all, active)
      msg <- glue(
        "Active GeneSetDb: {active_sets} of {total_sets} gene sets, ",
        "from {active_collections} / {total_collections} collections, ",
        "with {active_genes} / {total_genes} genes",
        active_sets = nrow(gs.active),
        total_sets = nrow(gs.all),
        active_collections = length(unique(gs.active[["collection"]])),
        total_collections = length(unique(gs.all[["collection"]])),
        active_genes = length(featureIds(gdb., active.only = TRUE)),
        total_genes = length(featureIds(gdb., active.only = FALSE)))
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
#' @importFrom multiGSEA.shiny userDefinedGeneSetDbUI
geneSetDbConfigUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    tags$div(id = ns("userdb"), userDefinedGeneSetDbUI(ns("user_gdb"))),
    tags$div(id = ns("msg"), uiOutput(ns("message"))))
}

#' @noRd
#' @export
initialized.ReactiveGeneSetDb <- function(x, ...) {
  is(x$gdb(), "GeneSetDb")
}

