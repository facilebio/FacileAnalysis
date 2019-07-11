# Interactivity and Vizualization over FacileDGEResults ========================

#' @noRd
#' @export
shine.FacileDgeAnalysisResult <- function(x, user = Sys.getenv("USER"),
                                          title = "Differential Expression Results",
                                          width = 800, height = 800,
                                          viewer = "pane", ...) {
  frunGadget(fdgeView, fdgeViewUI, x, dgeres = x, title = title,
             width = width, height = height, ...)
}

# @importFrom shiny callModule dialogViewer observeEvent runGadget stopApp
# @importFrom miniUI gadgetTitleBar miniContentPanel miniPage
# @importFrom FacileShine reactiveFacileDataStore
# shine.FacileDgeAnalysisResult <- function(x, user = Sys.getenv("USER"),
#                                           title = "Differential Expression Results",
#                                           width = 800, height = 800,
#                                           viewer = "pane", ...) {
#   ui <- miniPage(
#     useShinyjs(),
#     useSweetAlert(),
#     gadgetTitleBar(class(x)[1L]),
#     miniContentPanel(fdgeViewUI("view")),
#     NULL)
#
#   viewer <- gadget_viewer(viewer, title, width, height, ...)
#
#   server <- function(input, output, session) {
#     rfds <- ReactiveFacileDataStore(fds(x), "ds", user = user,
#                                     samples = samples(x))
#     # rfds <- callModule(reactiveFacileDataStore, "ds", fds(x), samples(x), user)
#     view <- callModule(fdgeView, "view", rfds, x, ...)
#     observeEvent(input$done, {
#       stopApp(invisible(NULL))
#     })
#     observeEvent(input$cancel, {
#       stopApp(invisible(NULL))
#     })
#   }
#
#   runGadget(ui, server, viewer = viewer, stopOnCancel = FALSE)
# }

#' @noRd
#' @export
#' @importFrom DT datatable formatRound
viz.FacileTtestAnalysisResult <- function(x, type = c("dge", "features"),
                                          ntop = 200, max_padj = 0.10,
                                          min_logFC = 1,
                                          features = NULL, round_digits = 3,
                                          event_source = "A", webgl = TRUE,
                                          ...) {
  type <- match.arg(type)
  treat_lfc <- x[["treat_lfc"]]
  if (!missing(min_logFC) && test_number(treat_lfc) && treat_lfc != min_logFC) {
    warning("DGE was run using TREAT. Minimum logFC is set to that threshold")
    min_logFC <- treat_lfc
  }

  fn <- if (type == "dge") .viz.dge_ttest else .viz.dge_features
  viz. <- fn(x, ntop = ntop, max_padj = max_padj, min_logFC = min_logFC,
             features = features, round_digits = round_digits,
             event_source = event_source, webgl = webgl, ...)

  if (type == "dge") {
    out <- bscols.(viz.[["title"]], viz.[["volcano"]], viz.[["datatable"]],
                   widths = c(12, 5, 7))
  } else {
    out <- viz.[["datatable"]]
  }

  out
}

#' @section Interacting with results:
#'
#' The `report` function will create an htmlwidget which can be explored by
#' the analyst or dropped into an Rmarkdown report.
#'
#' `report(result, "dge", max_padj = 0.05, min_logFC = 1)` will create a
#' side-by-side volcano and datatable for differential expression results.
#'
#' @export
#' @importFrom crosstalk bscols
#' @importFrom htmltools browsable tagList tags
#' @rdname fdge
report.FacileTtestAnalysisResult <- function(x, type = c("dge", "features"),
                                             ntop = 200, max_padj = 0.10,
                                             min_logFC = 1,
                                             features = NULL, round_digits = 3,
                                             event_source = "A", webgl = TRUE,
                                             caption = NULL, ...) {
  type <- match.arg(type)
  treat_lfc <- x[["treat_lfc"]]
  if (!missing(min_logFC) && test_number(treat_lfc) && treat_lfc != min_logFC) {
    warning("DGE was run using TREAT. Minimum logFC is set to that threshold")
    min_logFC <- treat_lfc
  }

  fn <- if (type == "dge") .viz.dge_ttest else .viz.dge_features
  viz. <- fn(x, ntop = ntop, max_padj = max_padj, min_logFC = min_logFC,
             features = features, round_digits = round_digits,
             event_source = event_source, treat_lfc = treat_lfc,
             webgl = webgl, ...)

  sdat <- viz.[["datatable"]][["data"]]
  mdef <- model(x)

  if (!is.null(caption)) {
    caption <- tags$p(caption)
  }

  if (type == "dge") {
    title <- viz.[["title"]]
    designf <- mdef[["design_formula"]]
    covariate <- param(mdef, "covariate")
    cstring <- mdef[["contrast_string"]]
    details <- tagList(
      tags$p(
        tags$strong("Design: "),
        tags$code(designf)),
      tags$p(
        tags$strong("Tested: "),
        tags$code(paste(covariate, cstring, sep = ": "))))

    header <- tagList(title, details, caption)
    out <- bscols.(header, viz.[["volcano"]], viz.[["datatable"]],
                   widths = c(12, 5, 7))
  } else {

  }

  out
}

# Internal =====================================================================

#' Extracts data and builds htmlwidgets to vizualize results of ttest
#'
#' @noRd
#' @importFrom crosstalk SharedData
#' @importFrom DT datatable formatRound
#' @importFrom htmltools tags
#' @importFrom plotly config layout plot_ly toWebGL
.viz.dge_ttest <- function(x, ntop = 200, max_padj = 0.10,
                           min_logFC = 1, round_digits = 3,
                           event_source = "A", treat_lfc = NULL,
                           webgl = TRUE, link_feature = TRUE,
                           server = FALSE, ...) {
  # Extract datatable
  if (!server && ntop > 1000) {
    warning("Non-serverside processing caps ntop to 1000 features",
            immediate. = TRUE)
    ntop <- ntop
  }

  dat.all <- result(x) %>%
    select(symbol, feature_id, logFC, padj, pval)

  dat.sig <- dat.all %>%
    filter(abs(logFC) >= min_logFC & padj <= max_padj)

  dat.up <- dat.sig %>%
    arrange(desc(logFC)) %>%
    head(ntop / 2)
  dat.down <- dat.sig %>%
    arrange(logFC) %>%
    head(ntop / 2)

  dat <- dat.up %>%
    bind_rows(dat.down) %>%
    rename(FDR = padj)

  # Add link for gene and axe feature_id
  # We assume ensembl ids for now
  is.ensid <- all(substr(dat[["feature_id"]], 1, 3) == "ENS")
  link_feature <- link_feature && is.ensid
  if (link_feature) {
    org <- gsub(" ", "_", organism(fds(x)))
    url <- sprintf("http://uswest.ensembl.org/%s/Gene/Summary?g=%s",
                   org, dat[["feature_id"]])
    html <- sprintf('<a href="%s" target="_blank">%s</a>',
                    url, dat[["symbol"]])
    dat <- mutate(dat, symbol = html, feature_id = NULL)
  }

  sdat <- SharedData$new(dat)

  # Generate volcano given the filtered data
  yaxis <- list(range = c(0, max(-log10(dat$pval))))
  p <- sdat %>%
    plot_ly(x = ~logFC, y = ~-log10(pval),
            type = "scatter", mode = "markers",
            hoverinfo = "text",
            # hoverlabel = list(bgcolor = "#3f3f3f"),
            hoverlabel = list(bgcolor = "#B3B2B3"),
            source = event_source,
            marker = list(color = "#838588"),
            text = ~paste0(
              "Symbol: ", symbol, "<br>",
              sprintf(paste0("logFC: %.", round_digits, "f<br>"), logFC),
              sprintf(paste0("FDR: %.", round_digits, "f<br>"), FDR),
              sprintf(paste0("pvalue: %.", round_digits, "f<br>"), pval))) %>%
    layout(yaxis = yaxis, dragmode = "select") %>%
    config(displaylogo = FALSE)

  if (webgl) p <- toWebGL(p)


  dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
  dtable <- sdat %>%
    datatable(filter = "top", extensions = "Scroller", style = "bootstrap",
              class = "compact", width = "100%", rownames = FALSE,
              options = dtopts, escape = !link_feature) %>%
    formatRound(c("logFC", "FDR", "pval"), round_digits)

  # Title
  .title <- "Differential Expression Results"
  if (!is.null(treat_lfc)) .title <- paste(.title, "(TREAT)")
  title <- tagList(
    tags$strong(.title),
    tags$br(),
    tags$span(sprintf("Top %d [FDR %.2f, abs(logFC) >= %.2f]",
                      nrow(dat), max_padj, min_logFC)))

  list(datatable = dtable, volcano = p, title = title)
}

.viz.dge_features <- function(x, ...) {

}
