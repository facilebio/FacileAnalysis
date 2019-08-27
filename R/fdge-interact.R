# Interactivity and Vizualization over FacileDGEResults ========================

#' For some reason, default views in a dialog or pane freeze up like y0!
#' @noRd
#' @export
shine.FacileDgeAnalysisResult <- function(x, user = Sys.getenv("USER"),
                                          title = "Differential Expression Results",
                                          viewer = "browser", ...) {
  frunGadget(fdgeView, fdgeViewUI, x, dgeres = x, title = title,
             viewer = viewer, ...)
}

#' The most common visualization downstream from a DGE analysis is either
#' a volcano plot, or the expression of signficant (or not) genes across the
#' samples tested in the DGE analysis.
#'
#' Therefore, viz(fdge) will return a volcano and viz(fdge, feature) will return
#' a gene-level expression (box) plot.
#'
#'
#' When we find interesting hits/genes downstream from a DGE analysis, I think
#' we most often want to query invidual hits to see why they came up signficant
#' (or not). Therefore the `vis.FcileDgeAnalysisResult` function will plot
#' gene-level expression across the groups tested in the `fdge` result,
#' for the gene(s) requested in the funciton call.
#'
#' @noRd
viz.FacileTtestAnalysisResult <- function(x, feature_id = NULL, type = NULL,
                                          highlight = NULL,
                                          round_digits = 3, event_source = "A",
                                          webgl = type %in% c("volcano", "ma"),
                                          ...) {
  try.tidy.class <- c("FacileFeatureSignature", "FacileFeatureRanks")
  if (test_multi_class(feature_id, try.tidy.class)) {
    feature_id <- tidy(feature_id)
  }
  if (is.data.frame(feature_id)) {
    feature_id <- feature_id[["feature_id"]]
  }

  if (is.null(type)) {
    if (is.null(feature_id) || length(feature_id) > 10L) {
      type <- "volcano"
    } else {
      type <- "feature"
    }
  }

  type <- match.arg(tolower(type), c("feature", "volcano", "ma"))

  if (type == "feature") {
    out <- .viz_dge_feature(x, feature_id, round_digits = round_digits,
                            event_source = event_source, ...)
  } else {
    force(webgl)
    out <- .viz_ttest_scatter(x, type, feature_id = feature_id,
                             highlight = highlight,
                             round_digits = round_digits,
                             event_source = event_source, webgl = webgl, ...)
  }

  out
}

#' @export
#' @noRd
viz.FacileAnovaAnalysisResult <- function(x, feature_id, round_digits = 3,
                                          event_source = "A", ...) {
  if (is.data.frame(feature_id)) {
    feature_id <- feature_id[["feature_id"]]
  }
  assert_string(feature_id)
  .viz_dge_feature(x, feature_id, round_digits = round_digits,
                   event_source = event_source, ...)
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

# Ttest Scatter Helper Functions ===============================================

#' Visualize the "global" scatter plot downstream of a DGE analysis, ie. the
#' volcano or MA plot.
.viz_ttest_scatter <- function(x, type = c("volcano", "ma"), feature_id = NULL,
                               ntop = 200,
                               round_digits = 3, event_source = "A",
                               width = NULL, height = NULL,
                               webgl = TRUE, dat = NULL,
                               color_aes = NULL,
                               hover = c("logFC", "pval", "padj", "symbol"),
                               highlight = NULL, ...) {
  type <- match.arg(type)
  if (is.null(dat)) {
    dat <- .data_ttest_scatter(x, feature_id = feature_id, ntop = ntop,
                               color_aes = NULL,
                               highlight = highlight, ...)
    color_aes <- attr(dat, "color_aes")
  }
  assert_subset(
    c("feature_id", "logFC", "pval", "AveExpr", "symbol"),
    colnames(dat))
  if (!is.null(color_aes)) assert_subset(color_aes, colnames(dat))
  if (is.character(hover)) {
    hover <- intersect(hover, colnames(dat))
  }

  if (type == "volcano") {
    dat <- mutate(dat, xval = logFC, yval = -log10(pval))
    xlabel <- "logFC"
    ylabel <- "-log10(pval)"
  } else {
    dat <- mutate(dat, xval = AveExpr, yval = logFC)
    xlabel <- "Average Expression"
    ylabel <- "logFC"
  }

  fp <- fscatterplot(
    dat, c("xval", "yval"), color_aes = color_aes,
    hover = hover, key = ".key", event_source = event_source,
    xlabel = xlabel, ylabel = ylabel,
    width = width, height = height, webgl = webgl, ...)
  fp
}

#' Extracts the data used for a "global view" of a DGE analysis, ie.
#' gene-level *statistics* for a comparison, which one might use in volcano or
#' MA plots.
#'
#' @noRd
.data_ttest_scatter <- function(x, feature_id = NULL, ntop = 200,
                                  color_aes = NULL, highlight = NULL,
                                  ...) {
  assert_class(x, "FacileTtestAnalysisResult")
  dat.all <- tidy(ranks(x))
  take.cols <- c("feature_id", "symbol", "name", "AveExpr",
                 "logFC", "pval", "padj", "t")
  dat <- dat.all[, intersect(take.cols, colnames(dat.all))]

  topn.fids <- c(
    head(dat.all[["feature_id"]], ntop / 2),
    tail(dat.all[["feature_id"]], ntop / 2))

  try.tidy.class <- c("FacileFeatureSignature", "FacileFeatureRanks")
  if (test_multi_class(highlight, try.tidy.class)) {
    highlight <- tidy(highlight)
  }
  if (is.data.frame(highlight)) {
    highlight <- highlight[["feature_id"]]
  }
  if (!is.character(highlight)) highlight <- NULL

  if (!is.character(feature_id)) {
    feature_id <- topn.fids
  }
  fids <- unique(c(feature_id, highlight))

  dat <- filter(dat, feature_id %in% fids)
  dat <- mutate(dat, .key = seq(nrow(dat)))

  if (is.character(highlight)) {
    dat[["highlight."]] <- ifelse(dat[["feature_id"]] %in% highlight,
                                  "highlight", "background")
    color_aes <- "highlight."
  }
  if (is.character(color_aes)) {
    assert_subset(color_aes, colnames(dat))
  }
  attr(dat, "color_aes") <- color_aes
  dat
}

# Feature Level Expression Helpers =============================================

#' Setup boxplot(s) for expression level data of gene(s) across the samples used
#' in a DGE analysis.
#'
#' @noRd
.viz_dge_feature <- function(x, feature_id, round_digits = 3, event_source = "A",
                             ylabel = NULL, title = NULL,
                             batch_correct = TRUE, prior.count = 0.1,
                             legendside = "bottom", width = NULL,
                             height = NULL, ...) {
  if (is.data.frame(feature_id)) {
    feature_id <- feature_id[["feature_id"]]
  }
  fid <- assert_character(feature_id, min.len = 1) # , max.len = 1)

  mod <- model(x)
  assay_name <- param(x, "assay_name")

  dat <- .data_dge_feature(x, fid, batch_correct = batch_correct,
                           prior.count = prior.count, ...)

  numer <- param(mod, "numer")
  denom <- param(mod, "denom")
  test.covariate <- param(mod, "covariate")
  batch <- param(mod, "batch")

  if (!test_string(ylabel)) {
    ylabel <- assay_units(fds(x), assay_name, normalized = TRUE)
    bc <- attr(dat, "batch_corrected")
    if (is.character(bc)) {
      ylabel <- paste(
        ylabel,
        sprintf("(batch corrected [%s])", paste(bc, collapse = ",")),
        sep = "<br>")
    }
  }
  if (!test_string(title) && !isFALSE(title)) {
    finfo <- features(fds(x), assay_name = assay_name, feature_id = fid)
    name <- finfo[["name"]]
    if (nchar(name) == 0 || is.na(name)) {
      title <- feature
    } else {
      title <- sprintf("%s (%s)", name, fid)
    }
  }
  if (isFALSE(title)) title <- NULL

  # TODO: Make this more elegant
  xaxis <- test.covariate
  xlabel <- test.covariate
  color.by <- test.covariate
  if (is.ttest(x)) {
    if (length(numer) + length(denom) > 2) {
      is.numer <- dat[[test.covariate]] %in% numer
      dat[["xaxe."]] <- ifelse(is.numer, "numer", "denom")
      xaxis <- "xaxe."
      xlabel <- "group"
    }
  }

  facet_aes <- if (length(fid) > 1) "feature_name" else NULL

  fplot <- fboxplot(dat, xaxis, "value", with_points = TRUE,
                    event_source = event_source, key = ".key",
                    color_aes = color.by, facet_aes = facet_aes,
                    hover = c("dataset", "sample_id", test.covariate, batch),
                    width = width, height = height, legendside = legendside,
                    xlabel = "", ylabel = ylabel, title = title)
  # TODO: Add DGE stats to the fdge feature plot. from dat:
  #   * [["pval"]]
  #   * [["padj"]]
  #   * [["logFC]] & ([["t"]])? OR [["F]]
  fplot
}

#' Extracts expresson data for a feature across the samples use for the test.
.data_dge_feature <- function(x, feature_id, batch_correct = TRUE,
                              prior.count = 0.1, ...) {
  # fid <- assert_string(feature_id)
  fid <- assert_character(feature_id, min.len = 1)

  mod <- model(x)
  test.covariate <- param(mod, "covariate")
  batch <- param(mod, "batch")
  corrected <- !is.null(batch) && isTRUE(batch_correct)

  if (!corrected) batch <- NULL
  assay_name <- param(x, "assay_name")
  samples. <- samples(x)

  dat <- fetch_assay_data(samples., fid, assay_name = assay_name,
                          normalized = TRUE, prior.count = prior.count,
                          batch = batch, main = test.covariate)

  dat <- droplevels(dat)
  dat[[".key"]] <- seq(nrow(dat))

  # transfer dge statistics to outgoing data.frame
  if (is.ttest(x)) {
    xfer.cols <- c("logFC", "pval", "padj", "t")
  } else {
    xfer.cols <- c("F", "pval", "padj")
  }

  feature <- filter(features(x), feature_id == fid)
  if (nrow(feature) == 0) {
    feature <- tibble()
  }

  stat.cols <- colnames(feature)
  for (xcol in xfer.cols) {
    if (xcol %in% stat.cols) {
      dat[[xcol]] <- feature[[xcol]]
    } else {
      dat[[xcol]] <- NA_real_
    }
  }

  attr(dat, "analysis_stats") <- intersect(xfer.cols, colnames(feature))
  attr(dat, "batch_corrected") <- if (corrected) batch else FALSE
  dat
}



#' @noRd
#' @export
#' @importFrom DT datatable formatRound
.viz.FacileTtestAnalysisResult <- function(x, type = c("dge", "features"),
                                           ntop = 200, max_padj = 0.10,
                                           min_logFC = 1,
                                           features = NULL, round_digits = 3,
                                           event_source = "A", webgl = TRUE,
                                           ...) {
  if (missing(type) && !is.null(features)) {
    # boxplot of features, otherwise
    type <- "features"
  } else {
    type <- match.arg(type)
  }
  if (type == "features") {
    stopifnot(!is.null(features))
    .viz.dge_features(x, features = features, ...)
  }

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

  dat.all <- tidy(x) %>%
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

