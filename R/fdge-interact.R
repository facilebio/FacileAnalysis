# Interactivity and Vizualization over FacileDGEResults ========================

#' Extracts data and builds htmlwidgets to vizualize results of ttest
#'
#' @noRd
#' @importFrom crosstalk SharedData
#' @importFrom DT datatable formatRound
#' @importFrom htmltools tags
.vizualize.dge_ttest <- function(x, ntop = 200, max_padj = 0.10,
                                 min_logFC = 1, round_digits = 3,
                                 event_source = "A", treat_lfc = NULL, ...) {
  # Extract datatable
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

  sdat <- SharedData$new(dat)

  # Generate volcano given the filtered data
  yaxis <- list(range = c(0, max(-log10(dat$pval))))
  p <- sdat %>%
    plot_ly(x = ~logFC, y = ~-log10(pval),
            type = "scatter", mode = "markers",
            hoverinfo = "text",
            text = ~paste0(
              "Symbol: ", symbol, "<br>",
              sprintf(paste0("logFC: %.", round_digits, "f<br>"), logFC),
              sprintf(paste0("FDR: %.", round_digits, "f<br>"), FDR),
              sprintf(paste0("pvalue: %.", round_digits, "f<br>"), pval))) %>%
    layout(yaxis = yaxis, dragmode = "select")

  dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
  dtable <- sdat %>%
    datatable(extensions = "Scroller", style = "bootstrap",
              class = "compact", width = "100%",
              options= dtopts) %>%
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

.vizualize.dge_features <- function(x, ...) {

}

#' @noRd
#' @export
#' @importFrom crosstalk bscols
#' @importFrom DT datatable formatRound
vizualize.FacileTtestDGEResult <- function(x, type = c("dge", "features"),
                                           ntop = 200, max_padj = 0.10,
                                           min_logFC = 1,
                                           features = NULL, round_digits = 3,
                                           event_source = "A", ...) {
  type <- match.arg(type)
  treat_lfc <- x[["treat_lfc"]]
  if (!missing(min_logFC) && test_number(treat_lfc) && treat_lfc != min_logFC) {
    warning("DGE was run using TREAT. Minimum logFC is set to that threshold")
    min_logFC <- treat_lfc
  }

  fn <- if (type == "dge") .vizualize.dge_ttest else .vizualize.dge_features
  viz <- fn(x, ntop = ntop, max_padj = max_padj, min_logFC = min_logFC,
            features = features, round_digits = round_digits,
            event_source = event_source, ...)

  if (type == "dge") {
    out <- suppressWarnings({
      bscols(viz[["title"]], viz[["volcano"]], viz[["datatable"]],
             widths = c(12, 4, 8))
    })
  } else {
    out <- viz[["datatable"]]
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
report.FacileTtestDGEResult <- function(x, type = c("dge", "features"),
                                        ntop = 200, max_padj = 0.10,
                                        min_logFC = 1,
                                        features = NULL, round_digits = 3,
                                        event_source = "A",
                                        caption = NULL, ...) {
  type <- match.arg(type)
  treat_lfc <- x[["treat_lfc"]]
  if (!missing(min_logFC) && test_number(treat_lfc) && treat_lfc != min_logFC) {
    warning("DGE was run using TREAT. Minimum logFC is set to that threshold")
    min_logFC <- treat_lfc
  }

  fn <- if (type == "dge") .vizualize.dge_ttest else .vizualize.dge_features
  viz <- fn(x, ntop = ntop, max_padj = max_padj, min_logFC = min_logFC,
            features = features, round_digits = round_digits,
            event_source = event_source, treat_lfc = treat_lfc, ...)

  sdat <- viz[["datatable"]][["data"]]
  mdef <- x[["model_def"]]

  if (type == "dge") {
    title <- viz[["title"]]
    details <- tagList(
      tags$p(
        tags$strong("Design: "),
        tags$code(mdef[["design_formula"]])),
      tags$p(
        tags$strong("Tested: "),
        tags$code(glue("{covariate}: ({numer}) - ({denom})", .envir = mdef))))

    header <- tagList(title, details, caption)
    out <- suppressWarnings({
      bscols(header, viz[["volcano"]], viz[["datatable"]],
             widths = c(12, 4, 8))
    })
  } else {

  }

  out
}
