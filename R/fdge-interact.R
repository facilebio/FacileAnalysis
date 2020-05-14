#' @include fdge.R
NULL

# Interactivity and Vizualization over FacileDGEResults ========================

#' For some reason, default views in a dialog or pane freeze up like y0!
#'
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
#' @export
viz.FacileTtestAnalysisResult <- function(x, features = NULL, type = NULL,
                                          highlight = NULL,
                                          round_digits = 3, event_source = "A",
                                          webgl = type %in% c("volcano", "ma"),
                                          ...) {
  try.tidy.class <- c("FacileFeatureSignature", "FacileFeatureRanks")
  if (test_multi_class(features, try.tidy.class)) {
    features <- tidy(features)
  }

  features <- extract_feature_id(features)
  highlight <- extract_feature_id(highlight)

  if (is.null(type)) {
    if (is.null(features) || length(features) > 10L) {
      type <- "volcano"
    } else {
      type <- "feature"
    }
  }

  type <- match.arg(tolower(type), c("feature", "volcano", "ma"))

  if (type == "feature") {
    out <- .viz_dge_feature(x, feature_id = features,
                            round_digits = round_digits,
                            event_source = event_source, ...)
  } else {
    force(webgl)
    out <- .viz_ttest_scatter(x, type, feature_id = features,
                             highlight_id = highlight,
                             round_digits = round_digits,
                             event_source = event_source, webgl = webgl, ...)
  }

  out
}

#' @export
#' @noRd
viz.FacileAnovaAnalysisResult <- function(x, features = NULL, round_digits = 3,
                                          event_source = "A", ...) {
  if (is.null(features)) {
    features <- slice(tidy(ranks(x)), 1L)
  }
  features <- extract_feature_id(features)
  .viz_dge_feature(x, features, round_digits = round_digits,
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
                                             features = NULL, highlight = NULL,
                                             round_digits = 3,
                                             event_source = "A", webgl = TRUE,
                                             caption = NULL, ...) {
  type <- match.arg(type)
  treat_lfc <- x[["treat_lfc"]]

  features.all <- features(x)

  if (nrow(features.all) <= ntop) {
    if (missing(max_padj)) max_padj <- 1
    if (missing(min_logFC)) min_logFC <- 0
  } else if (!missing(min_logFC) &&
             test_number(treat_lfc) &&
             treat_lfc != min_logFC) {
    warning("DGE was run using TREAT. Minimum logFC is set to that threshold")
    min_logFC <- treat_lfc
  }

  features <- extract_feature_id(features)
  highlight <- extract_feature_id(highlight)

  fn <- if (type == "dge") .viz.dge_ttest else .viz.dge_features
  viz. <- fn(x, ntop = ntop, max_padj = max_padj, min_logFC = min_logFC,
             features = features, highlight = highlight,
             round_digits = round_digits, event_source = event_source,
             treat_lfc = treat_lfc, webgl = webgl, ...)

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
#'
#' @noRd
#' @param feature_id we assume the parent exported functions already handled
#'   a features data.frame and passed in a character feature_id
.viz_ttest_scatter <- function(x, type = c("volcano", "ma"),
                               feature_id = NULL, ntop = 200,
                               highlight_id = NULL,
                               round_digits = 3, event_source = "A",
                               width = NULL, height = NULL,
                               webgl = TRUE, dat = NULL,
                               color_aes = NULL,
                               hover = c("logFC", "pval", "padj", "symbol")
                               , ...) {
  type <- match.arg(type)
  if (is.null(dat)) {
    dat <- .data_ttest_scatter(x, feature_id = feature_id, ntop = ntop,
                               highlight_id = highlight_id, color_aes = NULL,
                               ...)
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
    xlabel <- "log2FC"
    ylabel <- "-log10(pval)"
  } else {
    dat <- mutate(dat, xval = AveExpr, yval = logFC)
    xlabel <- "Average Expression"
    ylabel <- "log2FC"
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
#' @param feature_id,highlight we assume the parent exported functions already
#'   handled a features data.frame and passed in a character feature_id
.data_ttest_scatter <- function(x, feature_id = NULL, ntop = 200,
                                highlight_id = NULL, color_aes = NULL, ...) {
  assert_class(x, "FacileTtestAnalysisResult")
  dat.all <- tidy(ranks(x))
  take.cols <- c("feature_id", "symbol", "name", "AveExpr",
                 "logFC", "pval", "padj", "t")
  dat <- dat.all[, intersect(take.cols, colnames(dat.all))]

  if (is.null(feature_id)) {
    feature_id <- c(
      head(dat.all[["feature_id"]], ntop / 2),
      tail(dat.all[["feature_id"]], ntop / 2))
  } else {
    assert_character(feature_id, min.len = 1)
  }

  highlight_id <- extract_feature_id(highlight_id)

  fids <- unique(c(feature_id, highlight_id))

  dat <- filter(dat, feature_id %in% fids)
  dat <- mutate(dat, .key = seq(nrow(dat)))

  if (is.character(highlight_id)) {
    dat[["highlight."]] <- ifelse(dat[["feature_id"]] %in% highlight_id,
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
                             height = NULL,
                             legendtitle = NULL,
                             showticklabels = FALSE, ...) {
  fid <- extract_feature_id(feature_id)
  stopifnot(length(fid) > 0L)

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
    finfo <- features(fds(x), assay_name = assay_name, feature_ids = fid)
    name <- finfo[["name"]]
    if (nchar(name) == 0 || is.na(name)) {
      title <- fid
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
  if (!isFALSE(legendtitle)) {
    if (!test_string(legendtitle)) {
      legendtitle <- paste(param(model(x), "covariate"), collapse = ",")
    }

    fplot$plot <- layout(
      fplot$plot,
      legend = list(
        y = if (showticklabels) -0.1 else 0.05,
        title = list(
          text = sprintf("<b>%s</b>", legendtitle),
          font = list(size = 16))))
  }

  if (legendside == "bottom") {
    legend.y <- if (showticklabels) -0.1 else 0.05
  } else {
    legend.y <- 0.66
  }

  fplot$plot <- layout(
    fplot$plot,
    xaxis = list(showticklabels = showticklabels),
    legend = list(
      font = list(size = 16),
      y = legend.y))

  fplot
}

#' Extracts expression data for a feature across the samples use for the test.
#' @noRd
.data_dge_feature <- function(x, feature_id, batch_correct = TRUE,
                              prior.count = 0.1, ...) {
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
  if (is(x, "FacileTtestAnalysisResult")) {
    dat <- semi_join(dat, samples(x, tested_only = TRUE),
                     by = c("dataset", "sample_id"))
  }
  dat <- droplevels(dat)
  dat[[".key"]] <- seq(nrow(dat))

  # transfer dge statistics to outgoing data.frame
  if (is.ttest(x)) {
    xfer.cols <- c("logFC", "pval", "padj", "t")
  } else {
    xfer.cols <- c("F", "pval", "padj")
  }

  feature <- filter(features(x), feature_id %in% fid)
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
    features <- extract_feature_id(features)
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
#' @param highlight result of extract_feautre_id from parent call
.viz.dge_ttest <- function(x, ntop = 200, highlight = NULL,
                           max_padj = 0.10,
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
    filter(
      (abs(logFC) >= min_logFC & padj <= max_padj) |
        feature_id %in% highlight)

  dat.up <- dat.sig %>%
    arrange(desc(logFC)) %>%
    head(ntop / 2)
  dat.down <- dat.sig %>%
    arrange(logFC) %>%
    head(ntop / 2)

  dat <- dat.up %>%
    bind_rows(dat.down) %>%
    rename(FDR = padj) %>%
    distinct(feature_id, .keep_all = TRUE)

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

  if (is.character(highlight)) {
    dat[["highlight."]] <- ifelse(dat[["feature_id"]] %in% highlight,
                                  "fg", "bg")
  }

  sdat <- SharedData$new(dat)

  # Generate volcano given the filtered data
  yaxis <- list(range = c(0, max(-log10(dat$pval)) + 0.1))
  xaxis <- list(range = (max(abs(dat$logFC)) + 0.2) * c(-1, 1))
  p <- sdat %>%
    plot_ly(x = ~logFC, y = ~-log10(pval),
            type = "scatter", mode = "markers",
            hoverinfo = "text",
            # hoverlabel = list(bgcolor = "#3f3f3f"),
            color = if (is.character(highlight)) {
              ~highlight.
            }else {
              I("darkgrey")
            },
            colors = FacileViz::create_color_map(c("bg", "fg")),
            hoverlabel = list(bgcolor = "#B3B2B3"),
            source = event_source,
            # marker = list(color = "#838588"),
            text = ~paste0(
              "Symbol: ", symbol, "<br>",
              sprintf(paste0("logFC: %.", round_digits, "f<br>"), logFC),
              sprintf(paste0("FDR: %.", round_digits, "f<br>"), FDR),
              sprintf(paste0("pvalue: %.", round_digits, "f<br>"), pval))) %>%
    layout(yaxis = yaxis, xaxis = xaxis, dragmode = "select",
           showlegend = FALSE) %>%
    config(displaylogo = FALSE)

  if (webgl) p <- toWebGL(p)

  dtopts <- list(deferRender = TRUE, scrollY = 300, scroller = TRUE)
  noshow <- c("highlight")
  if (isTRUE(all.equal(dat$symbol, dat$feature_id))) {
    noshow <- c(noshow, "symbol")
  }
  dtable <- sdat %>%
    datatable(filter = "top", extensions = "Scroller", style = "bootstrap",
              class = "compact", width = "100%", rownames = FALSE,
              options = dtopts, escape = !link_feature,
              colnames = setdiff(colnames(dat), noshow)) %>%
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

