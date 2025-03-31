# Interactivity and vizualization over FacilePcaAnalysisResult =================

#' @noRd
#' @export
#' @examples
#' xpca <- fpca(FacileData::an_fds())
#' iviz <- viz(xpca, color_aes = "cell_type")
#' gg <- viz(xpca, color_aes = "cell_type", interactive = FALSE)
viz.FacilePcaAnalysisResult <- function(x, dims = NULL,
                                        type = c("scatter", "scree"),
                                        ..., 
                                        interactive = TRUE,
                                        color_aes = NULL,
                                        height = 400,
                                        width = 700,
                                        xlabel = "default",
                                        ylabel = "default",
                                        zlabel = "default",
                                        event_source = "A",
                                        webgl = FALSE) {
  type <- match.arg(type)
  if (type == "scatter") {
    res <- .viz_pca_scatter(x, dims = dims, color_aes = color_aes,
                            height = height, width = width, xlabel = xlabel,
                            ylabel = ylabel, zlabel = zlabel,
                            event_source = event_source, webgl = webgl, ...,
                            interactive = interactive)
  } else {
    res <- .viz_pca_scree(x, dims = dims, color_aes = color_aes,
                          height = height, width = width, xlabel = xlabel,
                          ylabel = ylabel, zlabel = zlabel,
                          event_source = event_source, webgl = webgl, ...,
                          interactive = interactive)
  }
  res
}

#' @noRd
#' @importFrom rlang .data
.viz_pca_scatter <- function(
    x, 
    dims = NULL,
    ..., 
    title = NULL,
    subtitle = NULL,
    color_aes = NULL, 
    color_map = NULL,
    shape_aes = NULL,
    shape_map = NULL,
    height = 400,
    width = 700, 
    xlabel = "default",
    ylabel = "default", 
    zlabel = "default",
    event_source = "A", 
    webgl = FALSE,
    interactive = TRUE,
    point_size = 2.5) {
  if (is.null(dims)) {
    dims <- 2
  }
  
  xx <- tidy(x)
  assert_integerish(dims, lower = 1L)
  if (length(dims) == 1L) {
    assert_int(dims, lower = 2, upper = 3L)
    dims <- 1:dims
  }
  assert_int(length(dims), lower = 1L, upper = 3L)
  dims <- unique(dims)

  pc.cols.all <- colnames(xx)[grep("^PC\\d+$", colnames(xx))]
  pc.cols.req <- paste0("PC", dims)
  pc.cols <- intersect(pc.cols.req, pc.cols.all)

  if (length(pc.cols) != length(dims)) {
    stop("There's something awry with the pc columns you requested")
  }

  # slimdown the data.frame to only include the PCs user asked for and rest
  # of the covariate data.
  xx.cols <- c(pc.cols, setdiff(colnames(xx), pc.cols.all))
  xx <- select(xx, {{xx.cols}})

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

  # TODO: Handle color by feature expression. Currently this is handled in the
  # fpcaView shiny module, but we want to push that down here so that users
  # using this function in an R session or for Rmd output can color by
  # expression as well.
  # https://github.com/facilebio/FacileAnalysis/issues/28
  if (!is.null(color_aes)) {
  }
  if (!is.null(shape_aes)) {
  }
  
  if (interactive) {
    p <- FacileViz::fscatterplot(
      xx, pc.cols, xlabel = xlabel, ylabel = ylabel,
      zlabel = zlabel, event_source = event_source,
      webgl = webgl, height = height, width = width,
      color_aes = color_aes, color_map = color_map,
      shape_aes = shape_aes, shape_map = shape_map, 
      ...)
  } else {
    x. <- pc.cols[[1]]
    y. <- pc.cols[[2]]
    
    p <- ggplot2::ggplot(xx) +
      ggplot2::aes(x = .data[[x.]], y = .data[[y.]]) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = xlabel,
        y = ylabel
      )
    if (!is.null(color_aes)) {
      checkmate::assert_choice(color_aes, colnames(xx))
      p <- p + ggplot2::aes(color = .data[[color_aes]])
      
      if (is.null(color_map)) {
        mapvals <- xx[[color_aes]]
        if (!is.numeric(mapvals)) {
          color_map <- FacileViz::create_color_map(xx[[color_aes]])
        }
      }
      p <- p + ggplot2::scale_color_manual(values = color_map)
    }
    if (!is.null(shape_aes)) {
      checkmate::assert_choice(shape_aes, colnames(xx))
      p <- p + ggplot2::aes(shape = .data[[shape_aes]])
      if (is.null(shape_map)) {
        mapvals <- xx[[shape_aes]]
        if (!is.numeric(mapvals)) {
          shape_map <- FacileViz::create_shape_map(mapvals)
          p <- p + ggplot2::scale_shape_manual(values = shape_map)
        }
      }
    }
  }
  p
}

#' @noRd
.viz_pca_scree <- function(x, dims = NULL, ..., title = "Scree Plot",
                           interactive = TRUE) {
  xx <- tidy(x)
  pc.cols.all <- colnames(xx)[grep("^PC\\d+$", colnames(xx))]
  pc.cols.req <- paste0("PC", dims)
  pc.cols <- intersect(pc.cols.req, pc.cols.all)
  
  xx.cols <- c(pc.cols, setdiff(colnames(xx), pc.cols.all))
  xx <- select(xx, {{xx.cols}})
  
  pcv <- x$percent_var * 100
  varexp <- dplyr::tibble(
    component = names(pcv),
    var_explained = unname(pcv),
    var_cumulative = cumsum(var_explained))
  
  fig <- plotly::plot_ly(data = varexp) |> 
    plotly::add_bars(x = ~component, y = ~var_explained, color = I("#1B62A6")) |> 
    plotly::add_lines(x = ~component, y = ~var_cumulative, color = I("#FD690F")) |> 
    plotly::add_markers(x = ~component, y = ~var_cumulative, color = I("#FD690F")) |> 
    plotly::layout(
      title = title,
      showlegend = FALSE,
      xaxis = list(
        title = "Component"),
      yaxis = list(
        title = "Variance explained (%)",
        range = list(0, 100)))
  
  out <- list(
    plot = fig,
    input_data = varexp,
    params = list())
  class(out) <- c("FacileScreePlot", "FacileViz")
  out
}

#' @export
#' @noRd
#' @importFrom shiny tagList tags
#' @param x The `FacilePcaAnalysisResult`
#' @param pcs The PC's to show (min 2, max 3). If a single integer is provided,
#'   PC's 1:`pcs` will be shown. If a vector is provided, then the PCs
#'   specified will be shown. Defaults is `3`, to show first 3 PCs
#' @param topn the number of top genes to enumerate that drive direction of PCs
#' @param feature_id_col the column name in `x[["row_covariates"]]` to use
#'   in the feature table.
report.FacilePcaAnalysisResult <- function(x, pcs = 3, with_features = TRUE,
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
#' @importFrom shiny tagList tags
.viz.fpca <- function(x, pcs = 3, ntop = 100, with_features = TRUE,
                      report_feature_as = NULL, webgl = FALSE, ...,
                      event_source = "A",
                      xlabel = "default",
                      ylabel = "default",
                      zlabel = "default") {
  scatter <- viz(x, pcs = pcs, ..., xlabel = xlabel, ylabel = ylabel,
                 zlabel = zlabel, event_source = event_source, webgl = webgl)
  # cnames <- colnames(scatter$input_data)
  cnames <- colnames(input_data(scatter)) # $input_data)
  pc.cols <- cnames[grepl("^PC\\d+", cnames)]

  pcv <- x$percent_var * 100

  if (with_features) {
    ranked <- ranks(x, type = "ranked", report_feature_as = report_feature_as)
    rtable <- ranked |>
      result() |>
      select(1, !!pc.cols) |>
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
