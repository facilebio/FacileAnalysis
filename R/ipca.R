#' Immersive PCA
#'
#' The code here is largely inspired by DESeq2's plotPCA.
#'
#' You should look at factominer:
#' * http://factominer.free.fr/factomethods/index.html
#' * http://factominer.free.fr/graphs/factoshiny.html
#'
#' @export
#' @rdname ipca
#'
#' @param x a data container
#' @return an ipca result
ipca <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("ipca")
}

#' @rdname ipca
#' @method ipca DGEList
#' @importFrom edgeR cpm
ipca.DGEList <- function(x, pcs = 1:10, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = TRUE, log = TRUE, ...) {
  m <- edgeR::cpm(x, prior.count = prior.count, log = log)
  ipca(m, pcs, ntop, row_covariates, col_covariates, ...)
}

#' @rdname ipca
#' @method ipca matrix
#' @importFrom matrixStats rowVars
ipca.matrix <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL,...) {
  assert_integerish(pcs, lower = 1L, upper = nrow(x))
  if (is(row_covariates, "data.frame")) {
    assert_true(nrow(x) == nrow(row_covariates))
    assert_character(rownames(row_covariates))
    assert_true(all(rownames(x) == rownames(row_covariates)))
  }
  if (is(col_covariates, "data.frame")) {
    assert_true(ncol(x) == nrow(col_covariates))
    assert_character(rownames(col_covariates))
    assert_true(all(colnames(x) == rownames(col_covariates)))
  }

  rv <- matrixStats::rowVars(x)
  take <- head(order(rv, decreasing = TRUE), ntop)

  pca <- prcomp(t(x[take,]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  names(percentVar) <- paste0("PC", seq(percentVar))

  dat <- as.data.frame(lapply(pcs, function(pc) pca$x[, pc]))
  colnames(dat) <- paste0("PC", pcs)
  rownames(dat) <- colnames(x)
  percentVar <- percentVar[colnames(dat)]

  if (is(col_covariates, "data.frame")) {
    dat <- cbind(dat, col_covariates[rownames(dat),,drop = FALSE])
  }

  result <- list(tidy = dat, percentVar = percentVar)
  class(result) <- c("ipca", "iresult")
  result
}

# iplot ========================================================================

#' @rdname ipca
#' @method iplot ipca
#' @export
iplot.ipca <- function(x, pcs = 1:3, color_aes = NULL, shape_aes = NULL,
                       hover = NULL, title = "Immersive PCA", ...) {
  xx <- tidy(x)
  assert_integerish(pcs, lower = 1L)
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

  if (is.character(color_aes)) {
    assert_subset(color_aes, colnames(xx))
    xx <- tidyr::unite_(xx, ".color", color_aes, remove = FALSE)
  } else {
    xx[[".color"]] <- I("black")
  }

  if (is.character(shape_aes)) {
    assert_subset(shape_aes, colnames(xx))
    xx <- tidyr::unite_(xx, ".shape", shape_aes, remove = FALSE)
  } else {
    xx[[".shape"]] <- I("circle")
  }

  if (is.character(hover)) {
    assert_subset(hover, colnames(xx))
    hvals <- lapply(hover, function(wut) {
      vals <- xx[[wut]]
      if (is.numeric(vals)) vals <- sprintf("%.2f", vals)
      if (!is.character(vals)) vals <- as.character(vals)
      paste0(wut, ": ", vals)
    })
    xx[[".hover"]] <- do.call(paste, c(hvals, list(sep = "<br>")))
  } else {
    xx[[".hover"]] <- ""
  }


  # recycle different shpaes to only sample from .shapes.all if there are
  # more unique levels in xx[[".shape"]] then we have shapes to choose from
  .shapes.all <- c("circle", "x", "o", "+")
  .shape.lvls <- unique(xx[[".shape"]])
  if (length(.shape.lvls) > length(.shapes.all)) {
    warning("There are more unique levels in shape aesthetic than shapes ",
            "to choose from. Shapes will be recycled", immediate. = TRUE)
  }
  .shapes.idx <- seq(length(.shape.lvls)) %% length(.shapes.all)
  .shapes.idx[.shapes.idx == 0] <- length(.shapes.all)
  .shapes <- .shapes.all[.shapes.idx]

  xf <- paste0("~", pc.cols[1])
  yf <- paste0("~", pc.cols[2])
  zf <- paste0("~", pc.cols[3])

  if (length(pcs) == 2L) {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], x$percentVar[pc.cols[1L]] * 100))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], x$percentVar[pc.cols[2L]] * 100))

    p <- plot_ly(xx, x = formula(xf), y = formula(yf), type = "scatter",
                 color = ~.color, mode = "markers",
                 symbol = ~.shape, symbols = .shapes,
                 text = ~.hover)
    p <- layout(p, xaxis = xaxis, yaxis = yaxis)
  } else {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], x$percentVar[pc.cols[1L]] * 100))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[3L], x$percentVar[pc.cols[3L]] * 100))
    zaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], x$percentVar[pc.cols[2L]] * 100))
    scene <- list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)

    p <- plot_ly(xx, x = formula(xf), z = formula(yf), y = formula(zf),
                 type = "scatter3d", mode = "markers",
                 color = ~.color,
                 symbol = ~.shape, symbols = .shapes,
                 text = ~.hover)
    p <- layout(p, scene = scene)
  }

  # TODO: Customize the hovertext for the plot using the `hover` argument
  p <- layout(p, title = title)
  p
}