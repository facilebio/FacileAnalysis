#' Reusable volcano plot
#'
#' Builds a Plotly volcano plot from a prepared differential-statistics table.
#' Data preparation and downstream interpretation stay with the caller.
#'
#' @param dat Data frame containing the x-axis, y-metric, and key columns.
#' @param ... Passed to [FacileViz::fscatterplot()].
#' @param x Column name for the x-axis value.
#' @param y Column name for the y-axis metric to transform with `-log10()`.
#' @param key Column name for the feature key.
#' @param labels Feature keys to label.
#' @param label_priority Ordered label columns to use before falling back to
#'   `key`.
#' @param hover Columns to include in Plotly hover text.
#' @param event_source Plotly event source.
#' @param label_limit Maximum number of labeled points.
#' @param x_cutoff,y_cutoff Optional cutoffs used to color points.
#' @param width,height,webgl Passed to [FacileViz::fscatterplot()].
#'
#' @return A Plotly htmlwidget.
#' @export
fvolcanoPlot <- function(
  dat,
  ...,
  x = "logFC",
  y = "pval",
  key = "feature_id",
  labels = NULL,
  label_priority = c("symbol", "feature_id"),
  hover = c("symbol", "logFC", "FDR"),
  event_source = NULL,
  label_limit = 20L,
  x_cutoff = NULL,
  y_cutoff = NULL,
  width = NULL,
  height = NULL,
  webgl = TRUE
) {
  assert_data_frame(dat)
  assert_string(x)
  assert_string(y)
  assert_string(key)
  assert_character(label_priority, min.len = 1L)
  assert_count(label_limit, positive = TRUE)
  .fvolcano_require_columns(dat, c(x, y, key))

  dat$yaxis <- .fvolcano_y_values(dat, y)
  label.dat <- .fvolcano_label_data(
    dat,
    dat$yaxis,
    labels,
    x,
    key,
    label_priority,
    label_limit
  )
  point.style <- .fvolcano_point_style(
    dat,
    dat$yaxis,
    x_cutoff,
    y_cutoff,
    labels,
    x,
    key,
    label_limit
  )
  axis.ranges <- .fvolcano_axis_ranges(dat, x)
  hover <- intersect(hover, colnames(dat))

  fplot <- FacileViz::fscatterplot(
    dat,
    c(x, "yaxis"),
    xlabel = x,
    ylabel = .fvolcano_y_label(y),
    width = width,
    height = height,
    hover = hover,
    webgl = webgl,
    event_source = event_source,
    key = key,
    ...
  )
  plt <- FacileViz::plot(fplot)
  plt <- .fvolcano_style_points(plt, point.style)
  plotly::layout(
    plt,
    annotations = .fvolcano_label_annotations(label.dat, dat, x),
    xaxis = list(range = axis.ranges$x),
    yaxis = list(range = axis.ranges$y)
  )
}

#' Volcano plot label choices
#'
#' @param dat Data frame containing feature labels and keys.
#' @param label_priority Ordered label columns to use before falling back to
#'   `key`.
#' @param key Column name for the feature key.
#'
#' @return A named character vector where names are display labels and values
#'   are feature keys.
#' @export
fvolcanoLabelChoices <- function(
  dat,
  label_priority = c("symbol", "feature_id"),
  key = "feature_id"
) {
  assert_data_frame(dat)
  assert_character(label_priority, min.len = 1L)
  assert_string(key)
  .fvolcano_require_columns(dat, key)
  .fvolcano_label_choices(dat, label_priority, key)
}

.fvolcano_y_choices <- function(y) {
  assert_character(y, min.len = 1L, any.missing = FALSE)
  if (is.null(names(y))) {
    names(y) <- y
  }
  y
}

.fvolcano_y_metric <- function(metric, y) {
  y <- .fvolcano_y_choices(y)
  if (length(metric) != 1L || !metric %in% unname(y)) {
    metric <- unname(y)[1L]
  }
  metric
}

.fvolcano_y_label <- function(metric) {
  sprintf("-log10(%s)", metric)
}

.fvolcano_y_values <- function(dat, metric) {
  -log10(dat[[metric]])
}

.fvolcano_feature_names <- function(dat, label_priority, key) {
  labels <- rep(NA_character_, nrow(dat))
  for (col in label_priority) {
    if (!col %in% colnames(dat)) {
      next
    }
    vals <- as.character(dat[[col]])
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- vals[missing]
  }
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- as.character(dat[[key]][missing])
  labels
}

.fvolcano_label_choices <- function(dat, label_priority, key) {
  labels <- .fvolcano_feature_names(dat, label_priority, key)
  duplicate.labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate.labels] <- paste0(
    labels[duplicate.labels],
    " (",
    dat[[key]][duplicate.labels],
    ")"
  )
  stats::setNames(as.character(dat[[key]]), labels)
}

.fvolcano_label_data <- function(
  dat, y.values, selected, x, key, label_priority, label_limit
) {
  if (is.null(selected)) {
    selected <- character()
  }
  selected <- utils::head(selected, label_limit)
  idx <- match(selected, as.character(dat[[key]]), nomatch = 0L)
  idx <- idx[idx > 0L]
  out <- dat[idx, , drop = FALSE]
  out$yaxis <- y.values[idx]
  out <- out[is.finite(out[[x]]) & is.finite(out$yaxis), , drop = FALSE]
  out$volcano_label <- .fvolcano_feature_names(out, label_priority, key)
  out
}

.fvolcano_axis_ranges <- function(dat, x) {
  padded_range <- function(vals, pad) {
    rng <- suppressWarnings(range(vals, finite = TRUE))
    if (!all(is.finite(rng))) {
      rng <- c(-0.5, 0.5)
    } else if (diff(rng) == 0) {
      rng <- rng + c(-0.5, 0.5)
    }
    span <- diff(rng)
    rng + c(-span * pad[1L], span * pad[2L])
  }
  list(
    x = padded_range(dat[[x]], c(0.08, 0.18)),
    y = padded_range(dat$yaxis, c(0.02, 0.14))
  )
}

.fvolcano_label_annotations <- function(label.dat, plot.dat, x) {
  if (nrow(label.dat) == 0L) {
    return(list())
  }

  lapply(seq_len(nrow(label.dat)), function(idx) {
    list(
      x = label.dat[[x]][idx],
      y = label.dat$yaxis[idx],
      xref = "x",
      yref = "y",
      text = htmltools::htmlEscape(label.dat$volcano_label[idx]),
      showarrow = TRUE,
      ax = 20,
      ay = -20,
      arrowhead = 0,
      arrowwidth = 1,
      arrowcolor = "#555555",
      bgcolor = "rgba(238, 238, 238, 0.86)",
      bordercolor = "#9a9a9a",
      borderwidth = 1,
      borderpad = 1,
      font = list(size = 10, color = "#222222"),
      opacity = 0.98,
      captureevents = FALSE
    )
  })
}

.fvolcano_point_palette <- list(
  background = "rgba(105, 105, 105, 0.18)",
  left = "rgba(0, 82, 204, 0.62)",
  right = "rgba(220, 35, 35, 0.62)",
  left_selected = "rgb(0, 82, 204)",
  right_selected = "rgb(220, 35, 35)",
  selected = "rgb(0, 0, 0)",
  outline = "rgba(35, 35, 35, 0.7)",
  no_outline = "rgba(0, 0, 0, 0)"
)

.fvolcano_point_style <- function(
  dat, y.values, x.cutoff, y.cutoff, selected = NULL, x, key, label_limit
) {
  colors <- rep(.fvolcano_point_palette$background, nrow(dat))
  line.colors <- rep(.fvolcano_point_palette$no_outline, nrow(dat))
  line.widths <- rep(0, nrow(dat))
  left.idx <- right.idx <- rep(FALSE, nrow(dat))
  x.cutoff <- suppressWarnings(as.numeric(x.cutoff))
  y.cutoff <- suppressWarnings(as.numeric(y.cutoff))
  has.cutoffs <- length(x.cutoff) == 1L && length(y.cutoff) == 1L &&
    is.finite(x.cutoff) && is.finite(y.cutoff)
  if (has.cutoffs) {
    x.cutoff <- abs(x.cutoff)
    left.idx <- dat[[x]] <= -x.cutoff & y.values >= y.cutoff
    right.idx <- dat[[x]] >= x.cutoff & y.values >= y.cutoff
    colors[left.idx] <- .fvolcano_point_palette$left
    colors[right.idx] <- .fvolcano_point_palette$right
  }
  selected.idx <- rep(FALSE, nrow(dat))
  if (length(selected)) {
    selected <- utils::head(selected, label_limit)
    selected.idx <- as.character(dat[[key]]) %in% selected
    colors[selected.idx & left.idx] <- .fvolcano_point_palette$left_selected
    colors[selected.idx & right.idx] <- .fvolcano_point_palette$right_selected
    colors[selected.idx & !left.idx & !right.idx] <-
      .fvolcano_point_palette$selected
  }
  outline.idx <- left.idx | right.idx | selected.idx
  line.colors[outline.idx] <- .fvolcano_point_palette$outline
  line.widths[outline.idx] <- 0.5
  list(color = colors, line.color = line.colors, line.width = line.widths)
}

.fvolcano_style_points <- function(plt, point.style) {
  if (length(plt$x$data) == 0L) {
    return(plt)
  }
  marker <- plt$x$data[[1L]]$marker
  if (is.null(marker)) {
    marker <- list()
  }
  marker$color <- point.style$color
  marker$size <- 4
  marker$opacity <- 1
  marker$line <- list(
    color = point.style$line.color,
    width = point.style$line.width
  )
  plt$x$data[[1L]]$marker <- marker
  plt
}

.fvolcano_require_columns <- function(dat, cols) {
  missing <- setdiff(cols, colnames(dat))
  if (length(missing)) {
    stop(
      sprintf(
        "Volcano data is missing required columns: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(dat)
}
