#' Immersive PCA
#'
#' The code here is largely inspired by DESeq2's plotPCA.
#'
#' You should look at factominer:
#' * http://factominer.free.fr/factomethods/index.html
#' * http://factominer.free.fr/graphs/factoshiny.html
#'
#' @export
#' @rdname fpca
#'
#' @param x a data container
#' @return an fpca result
fpca <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("fpca")
}

#' @rdname fpca
#' @method fpca DGEList
#' @importFrom edgeR cpm
fpca.DGEList <- function(x, pcs = 1:10, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = 3, log = TRUE, ...) {
  m <- edgeR::cpm(x, prior.count = prior.count, log = log)
  fpca(m, pcs, ntop, row_covariates, col_covariates, ...)
}

#' @rdname fpca
#' @method fpca matrix
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL, ...) {
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
  class(result) <- c("FacilePCA", "FacileAnalysisResult")
  result
}

# fplot ========================================================================

#' @rdname render
#' @method render FacilePCA
#'
#' @export
render.FacilePCA <- function(x, pcs = 1:3,
                             color_aes = NULL, color_map = NULL,
                             shape_aes = NULL, shape_map = NULL,
                             size_aes = NULL, size_map = NULL,
                             hover_aes = NULL, hover_map = NULL,
                             hover = NULL, ...) {
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

  p <- fscatterplot(xx, pc.cols,
                    color_aes = color_aes, color_map = color_map,
                    shape_aes = shape_aes, shape_map = shape_map,
                    size_aes = size_aes, size_map = size_map,
                    hover_aes = hover_aes, hover_map = hover_map,
                    hover = hover, ...)

  p$facile_analysis <- x

  if (length(pcs) == 2L) {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], x$percentVar[pc.cols[1L]] * 100))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], x$percentVar[pc.cols[2L]] * 100))
    p$plot <- layout(plot(p), xaxis = xaxis, yaxis = yaxis)
  } else {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], x$percentVar[pc.cols[1L]] * 100))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[3L], x$percentVar[pc.cols[3L]] * 100))
    zaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], x$percentVar[pc.cols[2L]] * 100))
    scene <- list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
    p$plot <- layout(plot(p), scene = scene)
  }

  class(p) <- c("FacilePCAPlot", class(p))
  p
}
