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
#'    [fdge_model_def()] function, and the shiny interface by the
#'    [fdgeModelDefRun()] module.
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
#' efds <- exampleFacileDataSet()
#' dge.crc <- efds %>%
#'   filter_samples(indication == "CRC") %>%
#'   fdgeGadget(viewer = "pane", with_volcano = FALSE)
#' dge.blca <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdgeGadget(viewer = "pane")
#' dge.comp <- compare(dge.crc, dge.blca)
#'
#' report(dge.comp)
#' shine(dge.comp)
#' }
fdgeGadget <- function(x, title = "Differential Expression Analysis",
                       height = 800, width = 1000, ...) {
  assert_multi_class(x, c("FacileDataStore", "facile_frame"))
  frunGadget(fdgeAnalysis, fdgeAnalysisUI, x, title = title,
             height = height, width = width, ...)
}

#' @rdname fdge-shiny
#' @section
#' Wrapper module to perform and interact with a differential expression result.
#'
#' This module can be embedded within a shiny app, or called from a gadget.
#'
#' @export
#' @importFrom shiny callModule
#' @importFrom shinyjs toggleElement
#' @param model_def the module returned from [fdgeModelDefModule()]
#' @return a `ShinyFacileDGEResult`, the output from [fdge()]
fdgeAnalysis <- function(input, output, session, rfds, with_volcano = TRUE,
                         ..., debug = FALSE) {
  model <- callModule(fdgeModelDefRun, "model", rfds, ..., debug = debug)
  dge <- callModule(fdgeRun, "dge", rfds, model, ..., debug = debug)
  view <- callModule(fdgeView, "view", rfds, dge, with_volcano = with_volcano,
                     ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = debug)

  # Only show the view UI when there is (at least) a defined model. When
  # done this way, the showSpinner() around the view's output datatable works
  # like a "wait, processing DGE" while it's running.
  observe({
    # model. <- model$result()
    # show <- is(model., "FacileDGEModelDefinition")
    res <- req(dge$result())
    show <- is(res, "FacileDGEResult")
    toggleElement("viewbox", condition = show)
  })

  vals <- list(
    result = dge,
    model = model,
    view = view,
    .ns = session$ns)
  class(vals) <- c("ShinyFdgeAnalysis", "ShinyFacileAnalysisResult")
  vals
}

#' @noRd
#' @export
result.ShinyFdgeAnalysis <- function(x, ...) {
  x[["result"]]$result()
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
fdgeAnalysisUI <- function(id, with_volcano = TRUE, ..., debug = FALSE,
                           bs4dash = isTRUE(getOption("facile.bs4dash"))) {
  ns <- NS(id)
  # box. <- if (bs4dash) bs4Dash::bs4Box else shinydashboard::box
  box. <- shinydashboard::box
  tagList(
    tags$div(id = ns("modelbox"),
             box.(title = "Model Definition", width = 12,
                  fdgeModelDefRunUI(ns("model"), debug = debug))),
    tags$div(id = ns("dgebox"),
             box.(title = "Testing Parameters", width = 12,
                  fdgeRunUI(ns("dge"), debug = debug))),
    hidden(tags$div(id = ns("viewbox"),
                    box.(title = NULL, #"Testing Results",
                         solidHeader = TRUE,
                         width = 12,
                         fdgeViewUI(ns("view"), with_volcano = with_volcano,
                                    debug = debug)))))
}

# Building block fdge shiny modules ============================================

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
#' @param model A linear model definition. Can be either a "naked"
#'   `FacileDGEModelDefinition` that is returned from a call to
#'   [fdge_model_def()], or the `ShinyDGEModelDefinition` object returned
#'   from the [fdgeModelDefRun()] module.
#' @param with_gsea Include option to run a GSEA?
#' @param ... passed into [fdge()]
#' @return A list of stuff. `$result` holds the `FacileDGEResult` wrapped in
#'   a `reactive()`, ie. a `ReactiveDGEResult`.
fdgeRun <- function(input, output, session, rfds, model, with_gsea = FALSE, ...,
                    debug = FALSE, .reactive = TRUE) {
  kosher <- is(model, "FacileDGEModelDefinition") ||
    is(model, "ShinyDGEModelDefinition")
  if (!kosher) {
    stop("Invalid object passed as `model`: ", class(model)[1L])
  }

  rmodel <- reactive({
    req(initialized(rfds))
    if (is(model, "FacileDGEModelDefinition")) {
      out <- model
    } else {
      out <- req(model$result())
    }
    out
  })

  assay <- callModule(assaySelect, "assay", rfds, .reactive = .reactive)

  # Sets up appropriate UI to retreive the params used to filter out
  # features measured using the selected assay
  ffilter <- callModule(fdgeFeatureFilter, "ffilter", rfds, model,
                        assay, ..., debug = debug)

  runnable <- reactive({
    model. <- rmodel()
    enable <- is(model., "FacileDGEModelDefinition") &&
      !is(model., "FacileFailedModelDefinition") &&
      !is(model., "IncompleteModelDefintion") &&
      initialized(assay)
  })

  # Update the dge_methods available given the selected assay
  observe({
    ainfo <- req(assay$assay_info())
    assay_type. <- ainfo$assay_type
    req(!unselected(assay_type.))
    method. <- input$dge_method
    methods. <- fdge_methods(assay_type.)$dge_method
    selected. <- if (method. %in% methods.) method. else methods.[1L]
    updateSelectInput(session, "dge_method", choices = methods.,
                      selected = selected.)
  })

  can_sample_weight <- reactive({
    method. <- input$dge_method
    req(!unselected(method.))
    fdge_methods() %>%
      filter(dge_method == method.) %>%
      pull(can_sample_weight)
  })
  sample_weights <- reactive({
    isTRUE(can_sample_weight() && input$sample_weights)
  })

  observe({
    toggleState("run", condition = runnable())
    toggleState("sample_weights", condition = can_sample_weight())
  })

  dge <- eventReactive(input$run, {
    req(runnable())
    model. <- rmodel()
    assay_name. <- assay$assay_info()$assay
    method. <- input$dge_method
    req(
      !unselected(assay_name.),
      !unselected(method.))
    sample_weights. <- sample_weights()
    treat_lfc. <- log2(input$treatfc)

    withProgress({
      fdge(model., assay_name = assay_name., method = method.,
           with_sample_weights = sample_weights., treat_lfc = treat_lfc.,
           filter = ffilter$method(),
           filter_min_count = ffilter$min_count(),
           filter_min_total_count = ffilter$min_total_count(),
           filter_min_expr = ffilter$min_expr())
    }, message = "Performing differential expression")
  })

  if (debug) {
    output$debug <- shiny::renderText({
      dge. <- req(dge())
      format(dge.)
    })
  }

  vals <- list(
    result = dge,
    .ns = session$ns)

  class(vals) <- "FacileShinyDGEResult"
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
#' @importFrom shinyWidgets dropdownButton prettyCheckbox
fdgeRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  out <- tagList(
    fluidRow(
      column(3, assaySelectUI(ns("assay"), choices = NULL)),
      column(
        1,
        tags$div(
          style = "padding-top: 1.7em",
          dropdownButton(
            inputId = ns("opts"),
            icon = icon("sliders"),
            status = "primary", circle = FALSE,
            width = "300px",
            selectInput(ns("dge_method"), label = "Method", choices = NULL),
            numericInput(ns("treatfc"), label = "Min. Fold Change",
                         value = 1, min = 1, max = 10, step = 0.25),
            prettyCheckbox(ns("sample_weights"), label = "Use Sample Weights",
                           status = "primary"),
            tags$h4("Filter Strategy"),
            wellPanel(
              fdgeFeatureFilterUI(ns("ffilter"), ..., debug = debug))))),
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


# Inner Helper Module for feature filtering parameters =========================

#' @noRd
#' @param assay_mod [FacileShine::assaySelect()] module, ie. an
#'   `AssaySelectInput` object.
fdgeFeatureFilter <- function(input, output, session, rfds, model,
                              assay_module, ..., debug = FALSE) {

  # TODO: Support returning a vector of feature id's
  filter_method <- reactive("default")

  min_count <- reactive({
    input$min_count
  })
  min_total_count <- reactive({
    input$min_total_count
  })
  min_expr <- reactive({
    input$min_expr
  })

  # Show/hide the appropropriate options. I would have used
  # shinyjs::toggleElement, but that's not working inside the dropdown
  # this is embedded in
  # observe({
  #   toggleElement("countopts", condition = count_options())
  #   toggleElement("expropts", condition = !count_options())
  # })
  output$show_count_options <- reactive({
    assay <- assay_module$assay_info()$assay_type
    if (assay %in% c("rnaseq", "umi", "isoseq")) "yes" else "no"
  })
  outputOptions(output, "show_count_options", suspendWhenHidden = FALSE)

  vals <- list(
    method = filter_method,
    min_count = min_count,
    min_total_count = min_total_count,
    min_expr = min_expr)
}

#' @noRd
#' @importFrom shiny conditionalPanel NS numericInput tagList tags
fdgeFeatureFilterUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)
  tagList(
    conditionalPanel(
      condition = "output.show_count_options == 'yes'", ns = ns,
      tags$div(id = ns("countopts"),
               numericInput(ns("min_count"), label = "Min. Count",
                            value = 10, min = 1, max = 50, step = 5),
               numericInput(ns("min_total_count"), label = "Min. Total Count",
                            value = 15, min = 1, max = 100, step = 10))),
    conditionalPanel(
      condition = "output.show_count_options == 'no'", ns = ns,
      tags$div(id = ns("expropts"),
               numericInput(ns("min_expr"), label = "Min. Expression",
                            value = 1, min = 0.25, max = 100, step = 5))))
}

# Visualize the fdge result ====================================================

#' Shiny module that interacts with result of fdgeRun module
#'
#' Provides an interactive view over the result of running the [fdgeRun()]
#' module.
#'
#' @export
#' @importFrom FacileViz fscatterplot fboxplot
#' @importFrom DT datatable dataTableProxy formatRound renderDT replaceData
#' @importFrom plotly event_data layout
#' @importFrom shiny downloadHandler req
#' @importFrom shinyjs toggleElement
#' @param result The result from calling [fdge()].
fdgeView <- function(input, output, session, rfds, result, with_volcano = TRUE,
                     ...,
                     feature_selection = session$ns("volcano"),
                     sample_selection = session$ns("samples"),
                     debug = FALSE) {

  # The FacileDGEResult object
  result. <- reactive({
    req(initialized(rfds))
    out <- if (is(result, "FacileDGEResult")) result else result$result()
    req(is(out, "FacileDGEResult"))
    out
  })

  assay_name. <- reactive(param(req(result.()), "assay_name"))
  samples. <- reactive(samples(req(result.()))) # should have covariate columns
  model. <- reactive(model(req(result.())))
  covariate. <- reactive(param(req(model.()), "covariate"))
  is.ttest <- reactive(is(req(result.()), "FacileTtestDGEResult"))

  observe({
    toggleElement("volcanobox", condition = is.ttest() && with_volcano)
    toggleElement("dlbuttonbox", condition = is(dge.stats.all(), "data.frame"))
  })

  # Two versions of the result table are stored:
  # 1. `dge.stats.all`: always the full table; and
  # 2. `dge.stats`: the subset of (1) which corresponds to the currently
  #    selected features in (1) that are brushed from the volcano. If no
  #    brushing is active, then (2) == (1)
  dge.stats.all <- reactive({
    .fdge <- req(result.())
    .result <- ranks(.fdge) %>% result()
    if (.fdge[["test_type"]] == "ttest") {
      out <- select(.result, symbol, feature_id, logFC, FDR = padj, pval)
    } else {
      out <- select(.result, symbol, feature_id, F, FDR = padj, pval)
    }
    out
  })

  # Visualization of DGE results -----------------------------------------------

  # Responds to selection events on the volcano plot. If no selection is active,
  # this is NULL, otherwise its the character vector of the "feature_id" values
  # that are brushed in the plot.
  selected_features <- reactive({
    selected <- event_data("plotly_selected", source = feature_selection)
    if (!is.null(selected)) selected <- selected$key
    selected
  })

  dge.stats <- reactive({
    stats. <- req(dge.stats.all())
    selected <- selected_features()
    if (!is.null(selected)) stats. <- filter(stats., feature_id %in% selected)
    stats.
  })

  if (with_volcano) {
    output$volcano <- renderPlotly({
      req(is.ttest())
      dat. <- req(dge.stats.all())
      dat.$yaxis <- -log10(dat.$pval)
      fplot <- fscatterplot(dat., c("logFC", "yaxis"),
                            xlabel = "logFC", ylabel = "-log10(pval)",
                            width = NULL, height = NULL,
                            hover = c("symbol", "logFC", "FDR"),
                            webgl = TRUE, event_source = feature_selection,
                            key = "feature_id")
      fplot$plot
    })
  }

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

  boxdata <- reactive({
    feature <- req(selected_result())
    dat <- req(samples.())
    scov <- covariate.()
    dat <- fetch_assay_data(rfds, feature, dat, normalized = TRUE,
                            prior.count = 0.25)
    if (!scov %in% colnames(dat)) {
      dat <- with_sample_covariates(dat, scov)
    }

    if (is.ttest()) {
      m <- req(model.())
      cov.vals <- c(param(m, "numer"), param(m, "denom"))
      keep <- dat[[scov]] %in% cov.vals # old school, !! isn't working
      dat <- dat[keep,,drop = FALSE]
    }

    dat <- droplevels(dat)
    mutate(dat, .key = seq(nrow(dat)))
  })

  output$boxplot <- renderPlotly({
    feature <- req(selected_result())
    dat <- req(boxdata())
    scov <- covariate.()
    if (is.ttest()) {
      m <- req(model.())
      numer <- param(m, "numer")
      denom <- param(m, "denom")
      if (length(numer) + length(denom) > 2) {
        is.numer <- dat[[scov]] %in% numer
        dat[["xaxe."]] <- ifelse(is.numer, "numer", "denom")
        xaxis <- "xaxe."
        xlabel <- "group"
        color.by <- scov
      } else {
        xaxis <- scov
        xlabel <- scov
        color.by <- scov
      }
    } else {
      xaxis <- scov
      xlabel <- scov
      color.by <- scov
    }

    fplot <- fboxplot(dat, xaxis, "value", with_points = TRUE,
                      event_source = sample_selection, key = ".key",
                      color_aes = color.by, hover = c(scov, "value"),
                      width = NULL, height = NULL, legendside = "bottom",
                      xlabel = "")
    out <- fplot$plot
    layout(out, title = sprintf("<b>%s</b>", feature[["symbol"]]))
  })

  # a dataset,sample_id tibble of the selected things, or NULL if none selected
  selected_samples <- reactive({
    selected <- event_data("plotly_selected", source = sample_selection)
    if (!is.null(selected)) {
      dat <- req(boxdata())
      selected <- select(filter(dat, .key == selected), dataset, sample_id)
    }
    selected
  })


  # Stats Table ----------------------------------------------------------------

  output$statsdl <- downloadHandler(
    filename = function() {
      res <- req(isolate(result.()))
      sprintf("DGE-%s_%s_%s.csv", param(res, "method"),
              res[["test_type"]], name(res))
    },
    content = function(file) {
      .fdge <- req(isolate(result.()))
      .result <- ranks(.fdge) %>% result()
      if ("padj" %in% colnames(.result)) {
        .result <- rename(.result, FDR = "padj")
      }
      c.order <- c("symbol", "name", "feature_id", "meta", "AveExpr",
                   "FDR", "pval")
      if (is(.fdge, "FacileAnovaDGEGadgetResult")) {
        c.order <- c(c.order, "F")
      } else {
        c.order <- c(c.order, "logFC", "CI.L", "CI.R", "B", "t")
      }
      c.order <- c(c.order, "assay")
      c.order <- intersect(c.order, colnames(.result))
      .result <- select(.result, !!c.order)
      write.csv(.result, file, row.names = FALSE)
    }
  )

  output$stats <- DT::renderDT({
    # Triggered when new result comes in
    res <- req(result.())
    dat <- req(dge.stats())
    num.cols <- colnames(dat)[sapply(dat, is.numeric)]

    dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
    dtopts <- list(pageLength = 15,
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
    selected_features = selected_features,
    selected_samples = selected_samples)
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
fdgeViewUI <- function(id, with_volcano = TRUE, ..., debug = FALSE) {
  ns <- NS(id)
  box <- shinydashboard::box

  fluidRow(
    column(
      width = 5,
      id = ns("vizcolumn"),
      # style = "padding: 2px; margin: 5px",
      # volcano
      if (with_volcano) {
        tags$div(id = ns("volcanobox"),
                 box(
                   width = 12,
                   withSpinner(plotlyOutput(ns("volcano")))))
      } else {
        NULL
      },
      # boxplot
      box(
        width = 12,
        id = ns("boxplotbox"),
        withSpinner(plotlyOutput(ns("boxplot"))))),
    column(
      width = 7,
      id = ns("statscolumn"),
      # style = "padding: 2px; margin: 5px",
      box(
        width = 12,
        id = ns("statbox"),
        title = "Differential Expression",
        tags$div(
          id = ns("dlbuttonbox"),
          style = "padding: 0 0 1em 0; float: right;",
          downloadButton(ns("statsdl"), "Download")),
        tags$div(style = "clear: right;"),
        withSpinner(DT::DTOutput(ns("stats"))))))
}
