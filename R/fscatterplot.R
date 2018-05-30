#' Draws interactive 2 or 3d scatterplots over data
#'
#' @rdname fscatterplot
#' @export
#' @param x a data object
#' @param axes the definition of the x,y,z axes
#' @return a plotly object
fscatterplot <- function(x, axes,
                         color_aes = NULL, color_map = NULL,
                         shape_aes = NULL, shape_map = NULL,
                         size_aes = NULL, size_map = NULL,
                         hover_aes = NULL, hover_map = NULL,
                         hover = NULL, ...) {
  UseMethod("fscatterplot", x)
}

fscatterplot.FacilePlot <- function(x, axes,
                                    color_aes = NULL, color_map = NULL,
                                    shape_aes = NULL, shape_map = NULL,
                                    size_aes = NULL, size_map = NULL,
                                    hover_aes = NULL, hover_map = NULL,
                                    hover = NULL, ...) {
  stop("This isn't real just yet, but it would look something like this")
  # delegate to fscatterpot.{tbl|data.frame}
  fscatterplot(tidy(x), axes,
               color_aes = color_aes, color_map = color_map,
               shape_aes = shape_aes, shape_map = shape_map,
               size_aes = size_aes, size_map = size_map,
               hover_aes = hover_aes, hover_map = hover_map,
               hover = hover, ...)
}

#' @rdname fscatterplot
#' @method fscatterplot data.frame
#' @export
fscatterplot.data.frame <- function(x, axes,
                                    color_aes = NULL, color_map = NULL,
                                    shape_aes = NULL, shape_map = NULL,
                                    size_aes = NULL, size_map = NULL,
                                    hover_aes = NULL, hover_map = NULL,
                                    face_aes = NULL, hover = NULL, ...) {
  assert_character(axes, min.len = 2L, max.len = 3L)
  assert_subset(axes, names(x))

  xx <- x
  xx <- with_color(xx, color_aes, aes_map = color_map, ...)
  xx <- with_shape(xx, shape_aes, aes_map = shape_map, ...)
  # xx <- with_size(xx, size_aes, ...)
  # xx <- with_hover(xx, hover_aes, ...)

  if (is.character(hover)) {
    assert_subset(hover, colnames(x))
    hvals <- lapply(hover, function(wut) {
      vals <- xx[[wut]]
      if (is.numeric(vals)) vals <- prettyNum(round(vals, 2), big.mark = ",")
      if (!is.character(vals)) vals <- as.character(vals)
      paste0(wut, ": ", vals)
    })
    xx[[".hover"]] <- do.call(paste, c(hvals, list(sep = "<br>")))
  } else {
    xx[[".hover"]] <- ""
  }

  xf <- paste0("~", axes[1])
  yf <- paste0("~", axes[2])
  zf <- paste0("~", axes[3])

  .colors <- suppressWarnings(aes_map(xx, "color"))
  if (is.null(.colors)) {
    .color <- I("black")
  } else {
    .color.columns <- attr(.colors, "columns")
    .color <- formula(paste0("~", .color.columns[["variable"]]))
  }

  .shapes <- suppressWarnings(aes_map(xx, "shape"))
  if (is.null(.shapes)) {
    .shape <- I("circle")
  } else {
    .shape.columns <- attr(.shapes, "columns")
    .shape <- formula(paste0("~", .shape.columns[["variable"]]))
  }

  if (length(axes) == 2L) {
    xaxis <- list(title = axes[1L])
    yaxis <- list(titl = axes[2L])
    p <- plot_ly(xx, x = formula(xf), y = formula(yf),
                 type = "scatter",  mode = "markers",
                 color = .color, colors = .colors,
                 symbol = .shape, symbols = .shapes,
                 text = ~.hover)
    p <- layout(p, xaxis = xaxis, yaxis = yaxis)
  } else {
    xaxis <- list(title = axes[1L])
    yaxis <- list(title = axes[3L])
    zaxis <- list(title = axes[2L])
    scene <- list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)

    p <- plot_ly(xx, x = formula(xf), z = formula(yf), y = formula(zf),
                 type = "scatter3d", mode = "markers",
                 color = .color, colors = .colors,
                 symbol = .shape, symbols = .shapes,
                 text = ~.hover,
                 marker = list(color = colors))
    p <- layout(p, scene = scene)
  }

  out <- list(plot = p, input_data = x, params = list())
  class(out) <- c("FacileScatterPlot", "FacilePlot")
  out
}

#' @rdname fscatterplot
#' @method fscatterplot tbl
#' @export
fscatterplot.tbl <- function(x, axes,
                             color_aes = NULL, color_map = NULL,
                             shape_aes = NULL, shape_map = NULL,
                             size_aes = NULL, size_map = NULL,
                             hover_aes = NULL, hover_map = NULL,
                             face_aes = NULL, hover = NULL, ...) {
  fscatterplot.data.frame(collect(x, Inf), axes,
                          color_aes = color_aes, color_map = color_map,
                          shape_aes = shape_aes, shape_map = shape_map,
                          size_aes = size_aes, size_map = size_map,
                          hover_aes = hover_aes, hover_map = hover_map,
                          hover = hover, ...)
}
