#' A ReactiveGeneSetDb for use in the shiiny world.
#'
#' This can be instantiated from a "static" or "reactive" GeneSetDb object.
#' It allows users to customize which genesets are active by:
#'
#' 1. Filtering out entire collections; and
#' 2. Filtering genesets based on min and max (gene) size.
#'
#' This should be in mutliGSEA.shiny, but is developed here for now to make
#' it easier to develop within the FacieAnalysis shiny world.
#'
#' @export
#' @importFrom multiGSEA geneSets
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shiny updateSliderInput
#'
#' @param gdb A static or reactive GeneSetDb object
#' @param min.gs.size,max.gs.size the default minimum and maximum geneset size
#'   set in the UI when `gdb` is first loaded or changes (when reactive)
reactiveGeneSetDb <- function(input, output, session, gdb,
                              min.gs.size = 2L, max.gs.size = Inf,
                              default_collections = NULL, ...) {
  assert_character(default_collections, min.len = 1L, null.ok = TRUE)
  assert_int(min.gs.size, lower = 1L, null.ok = TRUE)
  assert_int(max.gs.size, lower = 2L, null.ok = TRUE)

  state <- reactiveValues(
    gdb = NULL,
    collections = NULL,
    min.gs.size = min.gs.size,
    max.gs.size = max.gs.size)

  if (is(gdb, "reactive")) {
    observe({
      gdb. <- gdb()
      req(is(gdb., "GeneSetDb"))
      state$gdb <- gdb.
    })
  } else {
    state$gdb <- gdb
  }

  rgdb <- reactive(state$gdb)

  obsereve({
    gdb. <- req(rgdb())
    req(is(gdb., "GeneSetDb"))

    # ........................................................ collection picker
    gsets <- geneSets(gdb.)
    min.n <- min(gsets$N)
    max.n <- max(gsets$N)
    colls.all <- unique(gsets$collection)
    colls.selected <- intersect(colls.all, NULL)

    if (length(colls.selected) == 0L) {
      colls.selected <- colls.all
    }
    state$collections <- colls.selected
    updatePickerInput(session, "collections", choices = new.cols,
                      selected = colls.selected)

    # .............................................................. size slider
    if (!is.null(min.gs.size)) {
      min.val <- max(min.gs.size, min.n)
    } else {
      min.val <- min.n
    }
    if (!is.null(max.gs.size)) {
      max.val <- min(max.gs.size, max.n)
    } else {
      max.val <- max.n
    }
    state$min.gs.size <- min.val
    state$max.gs.size <- max.val
    updateSliderInput(session, "size", min = min.n, max = max.n,
                      value = c(min.n, max.n))
  })

  genesets <- reactive({
    gdb. <- req(rgdb())
    req(is(gdb., "GeneSetDb"))
    gsets <- geneSets(gdb.)
    if (is.character(state$collections)) {
      gsets <- filter(gsets, collection %in% stats$collections)
    }
    filter(gsets, N >= state$min.gs.size, N <= state$max.gs.size)
  })

  collections <- reactive({
    gs <- req(genesets())
    unique(gs$collection)
  })

  vals <- list(
    gdb = rgdb,
    geneSets = genesets,
    collections = collections,
    min.gs.size = reactive(state$min.gs.size),
    max.gs.size = reactive(state$max.gs.size),
    .state = state,
    .ns = session$ns)
  class(vals) <- c("ReactiveGeneSetDb")
  vals
}

#' @noRd
#' @importFrom shiny NS tagList sliderInput
#' @importFrom shinyWidgets pickerInput
reactiveGeneSetDbFilterUI <- function(id, min = 2, max = 100L, ...) {
  ns <- NS(id)

  tagList(
    pickerInput(ns("collections"), "Collections", choices = NULL,
                multiple = TRUE),
    sliderInput(ns("size"), "Set Size", min = 2, max = 100,
                value = c(min, max)),
  )
}
