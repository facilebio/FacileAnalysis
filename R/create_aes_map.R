# Helper functions to create categorical maps for shapes and colors

#' Maps colors to categorical values
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

#' Maps shaps to categorical values
create_shape_map <- function(vals, map = NULL) {
  stopifnot(is.categorical(vals))
  if (is.null(map)) {
    map <- 15:18
    map <- c(map, setdiff(1:25, map))
  }

  out <- xref.discrete.map.to.vals(map, vals)
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
#' @param map named character vector, where names are the entries found in
#'   `vals`
#' @param vals a categorical vector (character or factor)
#' @return a character vector like `map` but with recycled entries if the number
#'   of `length(unique(vals)) > length(map)`
xref.discrete.map.to.vals <- function(map, vals) {
  stopifnot(is.categorical(vals))
  stopifnot(is.character(map) || is.integerish(map))
  map.type <- if (is.character(map)) "char" else "int"

  if (is.factor(vals)) {
    uvals <- levels(vals)
  } else {
    uvals <- sort(unique(as.character(vals)))
  }

  if (is.null(names(map))) {
    out.map <- if (map.type == "char") character() else integer()
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
