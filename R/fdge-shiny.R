#' @include fdge.R
NULL

#' Fully interactive differential expression analysis
#'
#' Assembles a shiny UI to define all of the bits required to perform a
#' differential expression analysis over a predefined set of samples.
#'
#' An interactive differential expression analysis is divided into three
#' steps, each of which provides its own shiny module and interface. The
#' minimal input to this analysis is a pre-defined subset of samples to act on.
#'
#' These steps are:
#'
#' 1. Model matrix definition. The functionality is provided by the
#'    [flm_def()] function, and the shiny interface by the
#'    [flmDefRun()] module.
#' 2. Differential expression analysis. The functionality is
#'    defined by the [fdge()] function, and the shiny interface by the
#'    [fdgeRun()] module.
#' 3. Results display. The interactive display of the results is provided
#'    by the [fdgeView()] module.
#'
#' Th
#' @rdname fdge-shiny
#' @export
#'
#' @section Interactive Gadget:
#'
#' @examples
#' if (interactive()) {
#' # run tumor vs normal comparisons vs each, then run compare() on the results
#' options(facile.log.level.fshine = "trace")
#' efds <- FacileData::exampleFacileDataSet()
#' dge.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fdgeGadget(viewer = "pane")
#' dge.blca <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdgeGadget(viewer = "pane")
#' dge.comp <- compare(dge.crc, dge.blca)
#'
#' \dontrun{
#' tfds <- FacileDataSet("~/workspace/data/FacileData/dockerlink/FacileTcgaDataSet")
#' tsamples <- filter_samples(tfds, indication == "BRCA")
#' tdge <- fdgeGadget(tsamples, viewer = "browser")
#' }
#' report(dge.comp)
#' shine(dge.comp)
#' }
fdgeGadget <- function(x, title = "Differential Expression Analysis",
                       height = 800, width = 1000, viewer = "browser", ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fdgeAnalysis, fdgeAnalysisUI, x, title = title,
             height = height, width = width, viewer = viewer, ...)
}

#' @rdname fdge-shiny
#' @section Analysis Module:
#' Wrapper module to perform and interact with a differential expression result.
#'
#' This module can be embedded within a shiny app, or called from a gadget.
#'
#' @export
#' @importFrom shiny callModule
#' @importFrom shinyjs toggleElement
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ReactiveFacileDgeAnalysisResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, ..., debug = FALSE) {
  model <- callModule(flmDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)
  view <- callModule(fdgeView, "view", rfds, dge,  ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = debug)

  # Only show the view UI when there is (at least) a defined model. When
  # done this way, the showSpinner() around the view's output datatable works
  # like a "wait, processing DGE" while it's running.
  observe({
    res. <- req(faro(dge))
    show <- is(res., "FacileDgeAnalysisResult")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    main = dge,
    model = model,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ReactiveFacileDgeAnalysisResultContainer",
                   "ReactiveFacileAnalysisResultContainer")
  vals
}


#' @noRd
#' @export
#' @importFrom shiny
#'   fluidRow
#'   NS
#'   tagList
#'   tags
#'   wellPanel
#' @importFrom shinydashboard box
#' @importFrom shinyjs hidden
fdgeAnalysisUI <- function(id, ..., debug = FALSE,
                           bs4dash = isTRUE(getOption("facile.bs4dash"))) {
  ns <- NS(id)
  # box. <- if (bs4dash) bs4Dash::bs4Box else shinydashboard::box
  box. <- shinydashboard::box
  tagList(
    tags$div(id = ns("modelbox"),
             box.(title = "Model Definition", width = 12,
                  flmDefRunUI(ns("model"), debug = debug))),
    tags$div(id = ns("dgebox"),
             box.(title = "Testing Parameters", width = 12,
                  fdgeRunUI(ns("dge"), debug = debug))),
    hidden(tags$div(id = ns("viewbox"),
                    box.(title = NULL, #"Testing Results",
                         solidHeader = TRUE,
                         width = 12,
                         fdgeViewUI(ns("view"), debug = debug)))))
}

# Run ==========================================================================

#' Shiny module to configure and run fdge over a predefined model.
#'
#' @export
#' @importFrom shiny
#'   callModule
#'   eventReactive
#'   reactive
#'   renderText
#'   withProgress
#' @importFrom shinyjs toggleState
#' @importFrom FacileShine
#'   assaySelect
#'   unselected
#' @param model A linear model definition. Can be either an "innert"
#'   `FacileLinearModelDefinition` that is returned from a call to
#'   [flm_def()], or the `ReactiveDGEModelDefinition` object returned
#'   from the [flmDefRun()] module.
#' @param with_gsea Include option to run a GSEA?
#' @param ... passed into [fdge()]
#' @return A list of stuff. `$result` holds the `FacileDgeAnalysisResult`
#'   wrapped in a `reactive()`.
fdgeRun <- function(input, output, session, rfds, model, with_gsea = FALSE, ...,
                    debug = FALSE, .reactive = TRUE) {
  assert_class(rfds, "ReactiveFacileDataStore")
  assert_class(model, "FacileLinearModelDefinition")

  assay <- callModule(assaySelect, "assay", rfds, .reactive = .reactive)
  runopts <- callModule(fdgeRunOptions, "runopts", rfds, model, assay, ...)

  runnable <- reactive({
    initialized(model) && initialized(assay) && initialized(runopts)
  })

  observe({
    toggleState("run", condition = runnable())
  })

  dge <- eventReactive(input$run, {
    req(runnable())
    assay_name. <- assay$assay_info()$assay
    withProgress({
      fdge(model, assay_name = assay_name.,
           method = runopts$dge_method(),
           with_sample_weights = runopts$sample_weights(),
           treat_lfc = runopts$treat_lfc(),
           filter = runopts$feature_filter$method(),
           filter_min_count = runopts$feature_filter$min_count(),
           filter_min_total_count = runopts$feature_filter$min_total_count(),
           filter_min_expr = runopts$feature_filter$min_expr())
    }, message = "Performing differential expression")
  })

  if (debug) {
    output$debug <- shiny::renderText({
      dge. <- req(dge())
      format(dge.)
    })
  }

  vals <- list(
    faro = dge,
    .ns = session$ns)

  # Since `model` can be a reactive version of a FacileLinearModelDefinition, I
  # don't know how to get the Ttest or ANOVA state of that model without
  # being in a reactive context and, therefore, appending it to the class
  # of the outgoing result.
  class(vals) <- c("ReactiveFacileDgeAnalysisResult",
                   "ReactiveFacileAnalysisResult",
                   "FacileDgeAnalysisResult")
  vals
}

#' @noRd
#' @export
#' @importFrom FacileShine
#'   assaySelectUI
#' @importFrom shiny
#'   actionButton
#'   checkboxInput
#'   column
#'   fluidRow
#'   icon
#'   NS
#'   numericInput
#'   selectInput
#'   tagList
#'   tags
#'   verbatimTextOutput
#'   wellPanel
fdgeRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), choices = NULL)),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          fdgeRunOptionsUI(ns("runopts"), width = "300px"))),
      column(1, actionButton(ns("run"), "Run", style = "margin-top: 1.7em"))
    ))

  if (debug) {
    out <- tagList(
      out,
      tags$hr(),
      verbatimTextOutput(ns("debug"), placeholder = TRUE))
  }

  out
}

#' Runs fdge over a model definition created from a shiny gadget.
#'
#' @noRd
#' @export
fdge.ReactiveFacileLinearModelDefinition <- function(x, assay_name = NULL,
                                                     method = NULL,
                                                     filter = "default",
                                                     with_sample_weights = FALSE,
                                                     treat_lfc = NULL, ...) {
  fdge(faro(x), assay_name = assay_name, method = method, filter = filter,
       with_sample_weights = with_sample_weights, treat_lfc = treat_lfc, ...)
}


# View =========================================================================

#' Shiny module that interacts with result of fdgeRun module
#'
#' Provides an interactive view over the result of running the [fdgeRun()]
#' module.
#'
#' @export
#' @importFrom FacileViz fscatterplot fboxplot input_data
#' @importFrom DT datatable dataTableProxy formatRound renderDT replaceData
#' @importFrom plotly event_data layout
#' @importFrom shiny downloadHandler req
#' @importFrom shinyjs hide show toggleElement
#' @importFrom shinyWidgets updateSwitchInput
#' @param dgeres The result from calling [fdge()], or [fdgeRun()].
fdgeView <- function(input, output, session, rfds, dgeres, ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = FALSE) {
  state <- reactiveValues(
    # store brushed values from volcano plot here so that it can be reset
    # externally. When a new DGE is run, I want to turn off the volcano, and
    # blow out the selection. If we rely on just the reactive nature of
    # event_data("plotly_selected") then the selection would be stuck when
    # user turns off the volcano option
    volcano_select = tibble(assay_name = character(), feature_id = character()),
    boxplot_select = tibble(dataset = character(), sample_id = character()))

  # The FacileDGEResult object
  dge <- reactive({
    req(initialized(dgeres))
    faro(dgeres)
  })

  # When a new fdge result is produced, we may want to reset a few things.
  # Trigger that UI restting/cleanup work here.
  observeEvent(dge(), {
    # Hide volcanobox completely if new dge is not a ttest
    toggleElement("volcanobox", condition = is.ttest(dge()))
    # New results should turn off the volcano toggleSwitch if it is a ttest
    updateSwitchInput(session, "volcanotoggle", value = FALSE)
  }, priority = 5) # upping priority so some withProgress things hide quick

  # When the volcano boxplot is turned off (toggle is set to FALSE), nuke
  # any active selection
  observeEvent(input$volcanotoggle, {
    if (input$volcanotoggle) {
      show("volcanoplotdiv")
    } else {
      state$volcano_select <- tibble(feature_id = character())
      hide("volcanoplotdiv")
    }
  }, priority = 5)

  assay_name. <- reactive(param(req(dge()), "assay_name"))
  samples. <- reactive(samples(req(dge())))
  model. <- reactive(model(req(dge())))
  covariate. <- reactive(param(req(model.()), "covariate"))

  batch_corrected <- reactive({
    mod <- req(model.())
    batch <- param(mod, "batch")
    !is.null(batch)
  })

  # Two versions of the result table are stored:
  # 1. `dge.stats.all`: always the full table; and
  # 2. `dge.stats`: the subset of (1) which corresponds to the currently
  #    selected features in (1) that are brushed from the volcano. If no
  #    brushing is active, then (2) == (1)
  dge.stats.all <- reactive({
    dge. <- req(dge())
    ranked <- result(ranks(dge.))
    if (is.ttest(dge.)) {
      out <- select(ranked, symbol, feature_id, logFC, FDR = padj, pval)
    } else {
      out <- select(ranked, symbol, feature_id, F, FDR = padj, pval)
    }
    out
  })

  observe({
    toggleElement("dlbuttonbox", condition = is(dge.stats.all(), "data.frame"))
  })

  # Visualization of DGE results -----------------------------------------------

  output$volcano <- renderPlotly({
    plt <- NULL
    show <- input$volcanotoggle
    if (is.ttest(dge()) && show) {
      dat. <- req(dge.stats.all())
      if (debug) {
        # reduce volcano size
        dat. <- bind_rows(head(dat., 100), tail(dat., 100))
        dat. <- distinct(dat, feature_id, .keep_all = TRUE)
      }
      dat.$yaxis <- -log10(dat.$pval)
      fplot <- fscatterplot(dat., c("logFC", "yaxis"),
                            xlabel = "logFC", ylabel = "-log10(pval)",
                            width = NULL, height = NULL,
                            hover = c("symbol", "logFC", "FDR"),
                            webgl = TRUE, event_source = feature_selection,
                            key = "feature_id")
      plt <- plot(fplot)
    }
    plt
  })

  # Responds to selection events on the volcano plot. If no selection is active,
  # this is NULL, otherwise its the character vector of the "feature_id" values
  # that are brushed in the plot.
  observeEvent(event_data("plotly_selected", source = feature_selection), {
    selected <- event_data("plotly_selected", source = feature_selection)
    if (is.null(selected)) {
      selected <- character()
    } else {
      selected <- selected$key
    }
    if (!setequal(selected, state$volcano_select$feature_id)) {
      state$volcano_select <- tibble(
        assay_name = rep(assay_name.(), length(selected)),
        feature_id = selected)
    }
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  dge.stats <- reactive({
    stats. <- req(dge.stats.all())
    selected <- state$volcano_select
    if (nrow(selected)) {
      stats. <- filter(stats., feature_id %in% selected$feature_id)
    }
    stats.
  })

  # Visualization of selected genes across samples -----------------------------

  # The element in the statistics table that is selected, this will be an
  # assay/feature_id 1-row tibble, or NULL
  selected_result <- reactive({
    out <- input$stats_rows_selected
    aname <- isolate(assay_name.())
    if (!is.null(out)) {
      dat <- req(dge.stats())
      out <- tibble(assay = aname, feature_id = dat[["feature_id"]][out],
                    symbol = dat[["symbol"]][out])
    }
    out
  })

  # Keeps track of whether or not the model and fdge result used a "batch" term.
  # and shows/hides the UI element to view batch-corrected expression data
  # accordingly
  observe({
    enable_toggle <- batch_corrected()
    if (!enable_toggle) {
      updateSwitchInput(session, "batch_correct", value = FALSE)
    }
    toggleElement("batch_correct_container", condition = enable_toggle)
  })

  featureviz <- eventReactive(list(selected_result(), input$batch_correct), {
    dge. <- req(dge())
    bc <- batch_corrected() && input$batch_correct
    feature <- req(selected_result())
    viz(dge., feature, event_source = sample_selection, batch_correct = bc)
  })

  output$boxplot <- renderPlotly({
    plt <- req(featureviz())
    plot(plt)
  })

  output$boxplotdl <- downloadHandler(
    filename = function() {
      feature <- req(selected_result())
      req(nrow(feature) == 1L)
      symbol <- feature$symbol
      name <- "expression"
      if (!is.na(symbol) && !unselected(symbol)) {
        name <- paste(name, symbol, sep = "_")
      }
      name <- paste(name, feature[["feature_id"]], sep = "_")
      if (batch_corrected() && input$batch_correct) {
        name <- paste(name, "batch-corrected", sep = "_")
      }
      paste0(name, ".csv")
    },
    content = function(file) {
      dat <- req(featureviz()) %>%
        FacileViz::input_data() %>%
        select(dataset, sample_id, feature_id, feature_name, value) %>%
        with_sample_covariates()
      write.csv(dat, file, row.names = FALSE)
    }
  )

  # Responds to sample-selection events on the boxplot. If no selection is
  # active, this is NULL, otherwise it's <dataset>__<sample_id>-pastd
  # compound key
  observeEvent(event_data("plotly_selected", source = sample_selection), {
    selected <- event_data("plotly_selected", source = sample_selection)
    if (is.null(selected)) {
      selected <- tibble(dataset = character(), sample_id = character())
    } else {
      selected <- featureviz() %>%
        input_data() %>%
        slice(as.integer(selected$key)) %>%
        select(dataset, sample_id)
    }
    current <- paste(state$boxplot_select$dataset,
                     state$boxplot_select$sample_id)
    snew <- paste(selected$dataset, selected$sample_id)
    if (!setequal(current, snew)) {
      state$boxplot_select <- selected
    }
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # Stats Table ----------------------------------------------------------------

  output$statsdl <- downloadHandler(
    filename = function() {
      res <- req(isolate(dge()))
      test_type <- if (is.ttest(res)) "ttest" else "anova"
      sprintf("DGE-%s_%s_%s.csv", param(res, "method"), test_type, name(res))
    },
    content = function(file) {
      .fdge <- req(isolate(dge()))
      .result <- ranks(.fdge) %>% result()
      if ("padj" %in% colnames(.result)) {
        .result <- rename(.result, FDR = "padj")
      }
      c.order <- c("symbol", "name", "feature_id", "meta", "AveExpr",
                   "FDR", "pval")
      if (is.ttest(.fdge)) {
        c.order <- c(c.order, "logFC", "CI.L", "CI.R", "B", "t")
      } else {
        c.order <- c(c.order, "F")
      }
      c.order <- c(c.order, "assay")
      c.order <- intersect(c.order, colnames(.result))
      .result <- select(.result, !!c.order)
      write.csv(.result, file, row.names = FALSE)
    }
  )

  output$stats <- DT::renderDT({
    # Triggered when new result comes in
    res <- req(dge())
    dat <- req(dge.stats())
    num.cols <- colnames(dat)[sapply(dat, is.numeric)]

    dtopts <- list(deferRender = TRUE, scrollY = 300,
                   # scroller = TRUE,
                   pageLength = 15,
                   lengthMenu = c(15, 30, 50))

    dtable <- dat %>%
      datatable(filter = "top",
                # extensions = "Scroller",
                style = "bootstrap",
                class = "display", width = "100%", rownames = FALSE,
                selection = list(mode = "single", selected = 1, target = "row"),
                options = dtopts) %>%
      formatRound(num.cols, 3)
    dtable
  }, server = TRUE)
  # dtproxy <- dataTableProxy("stats")
  # observe({
  #   stats. <- req(dge.stats())
  #   replaceData(dtproxy, stats., resetPaging = FALSE)
  # })

  vals <- list(
    selected_features = reactive(state$volcano_select),
    selected_samples = reactive(state$boxplot_select),
    .ns = session$ns)
  return(vals)
}

#' @noRd
#' @export
#' @importFrom DT DTOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shiny column downloadButton fluidRow NS
#' @importFrom shinyjs hidden
#' @importFrom shinydashboard box
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets switchInput
fdgeViewUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  box <- shinydashboard::box

  volcano.box <- tags$div(
    style = "margin-top: 10px",
    id = ns("volcanobox"),
    box(
      width = 12, style = "padding-top: 10px",
      switchInput(ns("volcanotoggle"), size = "mini",
                  onLabel = "Yes", offLabel = "No", value = FALSE,
                  inline = TRUE),
      tags$span("Volcano Plot", style = "font-weight: bold"),
      tags$div(
        id = ns("volcanoplotdiv"),
        withSpinner(plotlyOutput(ns("volcano"))))))

  boxplot.box <- box(
    width = 12,
    id = ns("boxplotbox"),
    shinyjs::hidden(
      tags$div(
        id = ns("batch_correct_container"),
        switchInput(ns("batch_correct"), size = "mini",
                    onLabel = "Yes", offLabel = "No", value = TRUE,
                    inline = TRUE),
        tags$span("Batch Correction", style = "font-weight: bold"))),
    withSpinner(plotlyOutput(ns("boxplot"))))

  fluidRow(
    # Plots Column
    column(
      width = 5,
      id = ns("vizcolumn"),
      boxplot.box,
      downloadButton(ns("boxplotdl"), "Download Data"),
      volcano.box),
    # Stats Table Column
    column(
      width = 7,
      id = ns("statscolumn"),
      # style = "padding: 2px; margin: 5px",
      box(
        width = 12,
        id = ns("statbox"),
        title = "Statistics",
        tags$div(
          id = ns("dlbuttonbox"),
          style = "padding: 0 0 1em 0; float: right;",
          downloadButton(ns("statsdl"), "Download")),
        tags$div(style = "clear: right;"),
        withSpinner(DT::DTOutput(ns("stats"))))))
}
