# Interactivity and Vizualization over FacileDGEResults ========================

#' @noRd
#' @export
#' @importFrom crosstalk bscols SharedData
#' @importFrom DT datatable
vizualize.FacileTtestDGEResult <- function(x, result = c("features", "gsea"),
                                           with_boxplot = TRUE, ntop = 200,
                                           max_padj = 0.10, min_abs_logFC = 1,
                                           feature = NULL, event_source = "A",
                                           ...) {
  if (FALSE) {
    x <- FacileData::exampleFacileDataSet() %>%
     filter_samples(indication == "BLCA") %>%
     fdge_model_def(covariate = "sample_type",
                    numer = "tumor",
                    denom = "normal") %>%
      fdge(method = "voom", gsea = NULL)
    result <- "features"
    with_boxplot <- TRUE
    ntop <- 100
    max_padj <- 0.10
    min_abs_logFC <- 1
  }

  dat.all <- dge_stats(x) %>%
    select(symbol, feature_id, logFC, padj, pval, significant)

  dat.sig <- dat.all %>%
    filter(abs(logFC) >= min_abs_logFC & padj <= max_padj)

  dat.up <- dat.sig %>%
    arrange(desc(logFC)) %>%
    head(ntop / 2)
  dat.down <- dat.sig %>%
    arrange(logFC) %>%
    head(ntop / 2)

  dat <- bind_rows(dat.up, dat.down) %>%
    mutate(pval)
  sdat <- SharedData$new(dat)

  p <- sdat %>%
    plot_ly(x = ~logFC, y = ~-log10(pval), color = ~significant,
            type = "scatter", mode = "markers",
            hoverinfo = "text",
            text = ~paste("Symbol:", symbol,
                          sprintf("<br>logFC: %.2f", logFC),
                          sprintf("<br>FDR: %.2f", padj),
                          sprintf("<br>pval: %.2f", pval))) %>%
    layout(
      yaxis = list(range = c(0, max(-log10(dat$pval)))),
      dragmode = "select") %>%
    toWebGL()

  dtable <- datatable(
    sdat,
    extensions = "Scroller",
    style = "bootstrap", class = "compact", width = "100%",
    options=list(deferRender=TRUE, scrollY=300, scroller=TRUE))

  bscols(p, dtable, widths = c(4, 8))

}

#' @section Interacting with results:
#'
#' The `report` function will create an htmlwidget which can be explored by
#' the analyst or dropped into an Rmarkdown report.
#'
#' `report(result, "dge", max_padj = 0.05, min_abs_logFC = 1)` will create a
#' side-by-side volcano and datatable for differential expression results.
#'
#' @export
#' @rdname fdge
report.FacileTtestDGEResult <- function(x, result = "dge",
                                        max_padj = 0.01, min_abs_logFC = 1,
                                        event_source = "A", ...) {
  # result <- match.arg(result, c("dge", multiGSEA::resultNames(result$gsea)))
  result <- match.arg(result, "dge")

}
