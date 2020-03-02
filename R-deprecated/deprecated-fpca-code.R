# We don't need fpca.* methods for Bioc containers. We now rely on
# FacileBiocDataStore doing the right thing
#' @section DGEList
#' By default (`assay_name = "counts"`), the PCA will be performed on
#' log2 normalized counts. However, you may find that you'd like to store
#' something like a batch-corrected version of the counts back into the
#' DGEList, in which case you can provide the name of the element where this
#' matrix is stored as the `assay_name` parameter.
#'
#' @noRd
#' @export
#' @importFrom edgeR cpm
fpca.DGEList <- function(x, assay_name = NULL,
                         dims = min(5, ncol(x) - 1L),
                         filter = "default", ntop = 1000,
                         row_covariates = x$genes, col_covariates = x$samples,
                         batch = NULL, main = NULL, prior.count = 3, ...) {
  if (is.null(assay_name)) assay_name <- "counts"

  if (assay_name == "counts") {
    m <- edgeR::cpm(x, prior.count = prior.count, log = TRUE)
  } else {
    m <- x[[assay_name]]
    assert_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))
  }

  out <- fpca(m, dims, filter, ntop, row_covariates, col_covariates, batch,
              main, ...)

  out[["params"]][["assay_name"]] <- assay_name

  if ("dataset" %in% colnames(col_covariates)) {
    out[["samples"]][["dataset"]] <- col_covariates[["dataset"]]
  }
  if ("sample_id" %in% colnames(col_covariates)) {
    out[["samples"]][["sample_id"]] <- col_covariates[["sample_id"]]
  }

  out
}

#' @noRd
#' @export
fpca.EList <- function(x, assay_name = NULL,
                       dims = min(5, ncol(x) - 1L),
                       filter = "default", ntop = 1000,
                       row_covariates = x$genes, col_covariates = x$targets,
                       batch = NULL, main = NULL, ...) {
  if (is.null(assay_name)) assay_name <- "E"
  assert_choice(assay_name,  names(x))

  m <- x[[assay_name]]

  out <- fpca(m, dims, filter, ntop, row_covariates, col_covariates, batch,
              main, ...)
  out[["params"]][["assay_name"]] <- assay_name

  if ("dataset" %in% colnames(x$targets)) {
    out[["samples"]][["dataset"]] <- x$targets[["dataset"]]
  }
  if ("sample_id" %in% colnames(x$targets)) {
    out[["samples"]][["sample_id"]] <- x$targets[["sample_id"]]
  }

  out
}

#' @noRd
#' @export
fpca.ExpressionSet <- function(x, assay_name = NULL,
                               dims = min(5, ncol(x) - 1L),
                               filter = "default", ntop = 1000,
                               row_covariates = NULL, col_covariates = NULL,
                               batch = NULL, main = NULL, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase package required")

  if (is.null(row_covariates)) {
    row_covariates <- ns$fData(x)
  }
  if (is.null(col_covariates)) {
    col_covariates <- ns$pData(x)
  }

  if (is.null(assay_name)) {
    assay_name <- ns$assayDataElementNames(x)[[1L]]
  }
  m <- ns$assayDataElement(x, assay_name)
  assert_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))

  out <- fpca(m, dims, filter, ntop, row_covariates, col_covariates, batch,
              main, ...)

  out[["params"]][["assay_name"]] <- assay_name

  if ("dataset" %in% colnames(col_covariates)) {
    out[["samples"]][["dataset"]] <- col_covariates[["dataset"]]
  }
  if ("sample_id" %in% colnames(col_covariates)) {
    out[["samples"]][["sample_id"]] <- col_covariates[["sample_id"]]
  }

  out
}

#' This should be able to work on things like DESeqTransform objects, as well.
#' @noRd
#' @export
fpca.SummarizedExperiment <- function(x, assay_name = NULL,
                                      dims = min(5, ncol(x) - 1L),
                                      filter = "default", ntop = 1000,
                                      row_covariates = NULL,
                                      col_covariates = NULL,  batch = NULL,
                                      main = NULL, ...) {
  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  if (is.null(row_covariates)) {
    row_covariates <- ns4$as.data.frame.DataTable(ns$rowData(x))
  }
  if (is.null(col_covariates)) {
    col_covariates <- ns4$as.data.frame.DataTable(ns$colData(x))
  }

  if (is.null(assay_name)) {
    m <- ns$assays(x)[[1L]]
  } else {
    m <- ns$assay(x, assay_name)
  }
  assert_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))

  out <- fpca(m, dims, filter, ntop, row_covariates, col_covariates, batch,
              main, ...)

  out[["params"]][["assay_name"]] <- assay_name

  if ("dataset" %in% colnames(col_covariates)) {
    out[["samples"]][["dataset"]] <- col_covariates[["dataset"]]
  }
  if ("sample_id" %in% colnames(col_covariates)) {
    out[["samples"]][["sample_id"]] <- col_covariates[["sample_id"]]
  }

  out
}
