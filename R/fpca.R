#' Runs a facile PCA
#'
#' The code here is largely inspired by DESeq2's plotPCA.
#'
#' You should look at factominer:
#' * http://factominer.free.fr/factomethods/index.html
#' * http://factominer.free.fr/graphs/factoshiny.html
#'
#' @export
#' @importFrom multiGSEA eigenWeightedMean
#' @rdname fpca
#'
#' @param x a data container
#' @return an fpca result
fpca <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("fpca")
}

#' @export
#' @rdname fpca
#' @method fpca DGEList
#' @importFrom edgeR cpm
fpca.DGEList <- function(x, pcs = 1:10, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = 3, log = TRUE, ...) {
  m <- edgeR::cpm(x, prior.count = prior.count, log = log)
  fpca(m, pcs, ntop, row_covariates, col_covariates, ...)
}

#' @export
#' @rdname fpca
#' @method fpca matrix
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL, ...) {
  pcs.given <- !missing(pcs)
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
  assert_integerish(pcs, lower = 1L, upper = nrow(x))

  rv <- matrixStats::rowVars(x)
  take <- head(order(rv, decreasing = TRUE), ntop)

  xx <- x[take,,drop = FALSE]
  pca <- prcomp(t(xx))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  names(percentVar) <- paste0("PC", seq(percentVar))

  pcs.take <- paste0("PC", pcs)
  pcs.miss <- setdiff(pcs.take, colnames(pca$x))
  if (length(pcs.miss) && pcs.given) {
    warning("The following PCs were not included in the decomposition:\n  ",
            paste(pcs.miss, collapse = ", "))
    pcs.take <- intersect(pcs.take, colnames(pca$x))
  }

  dat <- as.data.frame(pca$x[, pcs.take])

  # Why was I doing this instead?
  # dat <- as.data.frame(lapply(pcs.take, function(pc) pca$x[, pc]))
  # colnames(dat) <- pcs.take
  # rownames(dat) <- colnames(x)

  percentVar <- percentVar[colnames(dat)]

  if (is(col_covariates, "data.frame")) {
    dat <- cbind(dat, col_covariates[rownames(dat),,drop = FALSE])
  }

  # Identify the percernt contribution each feature has to the PCs. It's running
  # a decomposition twice, but ...
  ewm <- eigenWeightedMean(xx, scale = FALSE)

  # Calculate correlation of each gene to PC1 -> PC4
  pc_cor <- sapply(paste0("PC", head(pcs, 4)), function(pc) {
    cor(t(xx), dat[[pc]])
  })

  rnames <- rownames(xx)
  if (is.null(rnames)) rnames <- as.character(seq(nrow(xx)))

  pc_cor <- bind_cols(
    tibble(row_name = rnames, row_idx = take),
    as.data.frame(pc_cor))

  result <- list(
    tidy = dat,
    factor_contrib = rename(ewm$factor.contrib, feature_id = "featureId"),
    percent_var = percentVar,
    pc_cor = pc_cor,
    taken = take)

  class(result) <- c("FacilePCAResult", "FacileReducedDimResult", "FacileAnalysisResult")
  result
}

#' @noRd
#' @export
#' @method print FacilePCAResult
print.FacilePCAResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacilePCAResult <- function(x, ...) {
  out <- paste(
    "===========================================================\n",
    sprintf("FacilePCAResult\n"),
    "-----------------------------------------------------------\n",
    "  n.observations\n",
    "  n.dimensions\n",
    "  % variance explained PC1 ... PCN\n",
    "===========================================================\n",
    sep = "")
  out
}
# fplot ========================================================================

#' @noRd
#' @method vizualize FacilePCAResult
#'
#' @export
vizualize.FacilePCAResult <- function(x, pcs = 1:3, ...) {
                                      # color_aes = NULL, color_map = NULL,
                                      # shape_aes = NULL, shape_map = NULL,
                                      # size_aes = NULL, size_map = NULL,
                                      # hover_aes = NULL, hover_map = NULL,
                                      # hover = NULL, ...) {
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

  # p <- fscatterplot(xx, pc.cols,
  #                   color_aes = color_aes, color_map = color_map,
  #                   shape_aes = shape_aes, shape_map = shape_map,
  #                   size_aes = size_aes, size_map = size_map,
  #                   hover_aes = hover_aes, hover_map = hover_map,
  #                   hover = hover, ...)

  p <- fscatterplot(xx, pc.cols, ...)

  p$facile_analysis <- x
  pcv <- x$percent_var * 100

  if (length(pcs) == 2L) {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], pcv[pc.cols[1L]]))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], pcv[pc.cols[2L]]))
    p$plot <- layout(plot(p), xaxis = xaxis, yaxis = yaxis)
  } else {
    xaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[1L], pcv[pc.cols[1L]]))
    yaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[3L], pcv[pc.cols[3L]]))
    zaxis <- list(title = sprintf("%s (%.2f%%)", pc.cols[2L], pcv[pc.cols[2L]]))
    scene <- list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
    p$plot <- layout(plot(p), scene = scene)
  }

  class(p) <- c("FacilePCAViz", "FacileReducedDimViz", class(p))
  p
}
