# Interactivity and vizualization over FacilePCAResults ========================

#' @noRd
#' @export
viz.FacilePCAResult <- function(x, pcs = 3, ...,
                                xlabel = "default",
                                ylabel = "default",
                                zlabel = "default",
                                event_source = "A",
                                webgl = FALSE) {
  xx <- result(x)
  assert_integerish(pcs, lower = 1L)
  if (length(pcs) == 1L) {
    assert_int(pcs, lower = 2, upper = 3L)
    pcs <- 1:pcs
  }
  assert_int(length(pcs), lower = 1L, upper = 3L)
  pcs <- unique(pcs)

  pc.cols.all <- colnames(xx)[grep("^PC\\d+$", colnames(xx))]
  pc.cols.req <- paste0("PC", pcs)
  pc.cols <- intersect(pc.cols.req, pc.cols.all)

  if (length(pc.cols) != length(pcs)) {
    stop("There's something awry with the pc columns you requested")
  }

  # slimdown the data.frame to only include the PCs user asked for and rest
  # of the covariate data.
  xx.cols <- c(pc.cols, setdiff(colnames(xx), pc.cols.all))
  xx <- xx[, xx.cols, drop = FALSE]

  pcv <- x$percent_var * 100

  if (xlabel == "default") {
    xlabel <- sprintf("%s (%.2f%%)", pc.cols[1L], pcv[pc.cols[1L]])
  }
  if (ylabel == "default") {
    ylabel <- sprintf("%s (%.2f%%)", pc.cols[2L], pcv[pc.cols[2L]])
  }
  if (zlabel == "default") {
    zlabel <- sprintf("%s (%.2f%%)", pc.cols[3L], pcv[pc.cols[3L]])
  }

  p <- fscatterplot(xx, pc.cols, xlabel = xlabel, ylabel = ylabel,
                    zlabel = zlabel, event_source = event_source,
                    webgl = webgl, ...)
  p
}

#' @export
#' @noRd
#' @importFrom shiny tagList tags
#' @param x The `FacilePCAResult`
#' @param pcs The PC's to show (min 2, max 3). If a single integer is provided,
#'   PC's 1:`pcs` will be shown. If a vector is provided, then the PCs
#'   specified will be shown. Defaults is `3`, to show first 3 PCs
#' @param topn the number of top genes to enumerate that drive direction of PCs
#' @param feature_id_col the column name in `x[["row_covariates"]]` to use
#'   in the feature table.
report.FacilePCAResult <- function(x, pcs = 3, with_features = TRUE,
                                   ntop = 100, report_feature_as = NULL,
                                   caption = NULL, ...,
                                   event_source = "A",
                                   xlabel = "default",
                                   ylabel = "default",
                                   zlabel = "default",
                                   webgl = FALSE) {
  viz. <- .viz.fpca(x, pcs = pcs, ntop = ntop,
                    with_features = with_features,
                    report_feature_as = report_feature_as, ...,
                    event_source = event_source, xlabel = xlabel,
                    ylabel = ylabel, zlabel = zlabel)

  if (!is.null(caption)) {
    caption <- tags$p(caption)
  }

  header <- tagList(viz.[["title"]], caption)

  if (with_features) {
    out <- bscols.(header, viz.[["plot"]]$plot, viz.[["datatable"]],
                   widths = c(12, 7, 5))
  } else {
    out <- bscols.(header, viz.[["plot"]], widths = c(12, 12))
  }

  out
}

# Internal =====================================================================

#' Internal worker vizualization function
#' @noRd
#' @importFrom DT datatable
#' @importFrom FacileViz fscatterplot
#' @importFrom shiny tagList tags
.viz.fpca <- function(x, pcs = 3, ntop = 100, with_features = TRUE,
                      report_feature_as = NULL, webgl = FALSE, ...,
                      event_source = "A",
                      xlabel = "default",
                      ylabel = "default",
                      zlabel = "default") {
  scatter <- viz(x, pcs = pcs, ..., xlabel = xlabel, ylabel = ylabel,
                 zlabel = zlabel, event_source = event_source, webgl = webgl)
  cnames <- colnames(scatter$input_data)
  pc.cols <- cnames[grepl("^PC\\d+", cnames)]

  pcv <- x$percent_var * 100

  if (with_features) {
    ranked <- ranks(x, type = "ranked", report_feature_as = report_feature_as)
    rtable <- ranked %>%
      select(1, !!pc.cols) %>%
      head(ntop)
    dtable <- datatable(
      rtable, extensions = "Scroller", style = "bootstrap",
      class = "compact", width = "100%", rownames = FALSE,
      options = list(deferRender = TRUE, scrollY = 300, scroller = TRUE))
  } else {
    dtable <- NULL
  }

  nfeatures <- nrow(x[["factor_contrib"]])

  title <- tagList(
    tags$strong("PCA Dimensionality Reduction"),
    tags$br(),
    tags$span(sprintf("Top %d features", nfeatures)),
    tags$br(),
    tags$p("Variance explained per PC:"),
    tags$ul(
      style = "list-style-type: none",
      lapply(pc.cols, function(pc) {
        tags$li(style = "float: left; margin-right: 1em",
                sprintf("%s: %.02f%%", pc, pcv[pc]),
                escape = FALSE)
      })
    ),
    tags$div(style = "clear: left"))
  out <- list(datatable = dtable, plot = scatter, title = title, pcv = pcv)
  out
}
