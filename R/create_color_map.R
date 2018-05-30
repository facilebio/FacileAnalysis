#' Maps colors to values
#'
#' @section Categorical Map:
#' Map can be a RColorBrewer name, or a vector of colors. Colors will be
#' recycled if `length(map) <`
#' If `vals`, map can be is categoricasl
create_color_map <- function(vals, map = NULL) {
  is.cat <- is.categorical(vals)

  if (is.cat) {
    if (is.null(map)) map <- mucho.colors()
    if (is.brewer.map.name(map)) {
      map <- suppressWarnings(brewer.pal(20, map))
    }
    if (!is.character(map)) {
      stop("The color map should be a list of colors by now")
    }
    out <- xref.discrete.map.to.vals(map, vals)
  } else {
    stop("Not mapping real values yet")
  }

  out
}

#' @noRd
#' @importFrom RColorBrewer brewer.pal.info
is.brewer.map.name <- function(x) {
  is.character(x) && length(x) == 1L && x %in% rownames(brewer.pal.info)
}

#' @noRd
#' @importFrom RColorBrewer brewer.pal
mucho.colors <- function() {
  s1 <- RColorBrewer::brewer.pal(9, "Set1")
  s2 <- RColorBrewer::brewer.pal(8, "Set2")
  s3 <- RColorBrewer::brewer.pal(12, "Set3")
  muchos <- c(s1, s2[1:8])
}

#' @noRd
xref.discrete.map.to.vals <- function(map, vals) {
  stopifnot(is.categorical(vals))
  stopifnot(is.character(map))

  if (is.factor(vals)) {
    uvals <- levels(vals)
  } else {
    uvals <- sort(unique(as.character(vals)))
  }

  if (is.null(names(map))) {
    out.map <- character()
    rest.map <- map
  } else {
    out.map <- map[names(map) %in% uvals]
    rest.map <- unname(map[!names(map) %in% names(out.map)])
  }


  remain <- setdiff(uvals, names(out.map))
  if (length(remain)) {
    cols <- unname(c(rest.map, out.map))
    idxs <- seq(remain) %% length(cols)
    idxs[idxs == 0] <- length(cols)
    rest.map <- cols[idxs]
    names(rest.map) <- remain
    out.map <- c(out.map, rest.map)
  }

  out.map
}
