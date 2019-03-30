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
#' @examples
#' # Generate a data matrix with a sister covariate table
#' mat <- local({
#'   m <- matrix(rnorm(100 * 10), nrow = 100)
#'   colnames(m) <- letters[1:10]
#'   rownames(m) <- head(unique(
#'     replicate(200, paste(sample(letters, 5), collapse = ""))),
#'     nrow(m))
#'   m
#' })
#'
#' # generate some sample (column) covariates
#' pdat <- local({
#'   pd <- FacileAnalysis:::example_aes_data_table(10, n.cats = 5)
#'   pd <- as.data.frame(pd)
#'   rownames(pd) <- colnames(mat)
#'   pd
#' })
#'
#' # analyze and vizualize
#' res <- fpca(mat, col_covariates = pdat)
#' viz(res, color_aes = "category")
fpca <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("fpca", x)
}

#' @export
#' @rdname fpca
#' @importFrom edgeR cpm
fpca.DGEList <- function(x, pcs = 1:10, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = 3, log = TRUE, ...) {
  m <- edgeR::cpm(x, prior.count = prior.count, log = log)
  fpca(m, pcs, ntop, row_covariates, col_covariates, ...)
}

#' @export
#' @rdname fpca
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL, ...) {
  pcs.given <- !missing(pcs)
  assert_integerish(pcs, lower = 1L, upper = nrow(x))

  if (is.null(rownames(x))) rownames(x) <- as.character(seq(nrow(x)))
  if (is.null(row_covariates)) {
    row_covariates <- data.frame(symbol = rownames(x), row.names = rownames(x),
                                 stringsAsFactors = FALSE)
  }
  assert_data_frame(row_covariates)
  assert_true(nrow(x) == nrow(row_covariates))
  assert_character(rownames(row_covariates))
  assert_true(all(rownames(x) == rownames(row_covariates)))

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

  # Identify the percernt contribution each feature has to the PCs.
  # We are the same decomposition twice, but convenience wins for now.
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
    result = dat,
    factor_contrib = rename(ewm$factor.contrib, feature_id = "featureId"),
    percent_var = percentVar,
    pc_cor = pc_cor,
    row_covariates = row_covariates,
    taken = take)

  class(result) <- c("FacilePCAResult", "FacileReducedDimResult",
                     "FacileAnalysisResult")
  result
}

#' Extracts the genes that contribute most to each PC
#'
#' @export
#' @noRd
#' @param report the column in x$row_covariates to use as the listing of the
#'   feature in the table of ranks
ranks.FacilePCAResult <- function(x, report = NULL, ...) {
  fcontrib <- x[["factor_contrib"]]
  rdata <- x[["row_covariates"]]

  if (is.null(report)) report <- colnames(rdata)[1L]
  assert_choice(report, colnames(rdata))

  pc.cols <- colnames(fcontrib)[grepl("PC\\d+", colnames(fcontrib))]
  out <- tibble(rank = seq(nrow(fcontrib)))
  for (pc in pc.cols) {
    o <- order(fcontrib[[pc]], decreasing = TRUE)
    xref <- match(fcontrib[["feature_id"]][o], rownames(rdata))
    out[[pc]] <- rdata[[report]][xref]
  }
  out
}

#' @noRd
#' @export
print.FacilePCAResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacilePCAResult <- function(x, ...) {
  n.features <- nrow(x[["factor_contrib"]])
  pcv <- x$percent_var * 100
  pcvs <- paste(names(pcv), sprintf("%.2f%%", pcv), sep = ": ")
  pcvu <- paste(head(pcvs, 5), collapse = "\n  ")

  out <- paste(
    "===========================================================\n",
    sprintf("FacilePCAResult\n"),
    "-----------------------------------------------------------\n",
    "Number of features used:", n.features, "\n",
    "Number of PCs:", length(pcv), "\n",
    "Variance explained:\n  ", pcvu, "\n",
    "===========================================================\n",
    sep = "")
  out
}

