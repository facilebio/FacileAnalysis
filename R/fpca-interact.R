# Interactivity and vizualization over FacilePCAResults ========================

#' @noRd
#' @method vizualize FacilePCAResult
#'
#' @export
#' @importFrom crosstalk bscols SharedData
#' @importFrom DT datatable
#'
#' @param x The `FacilePCAResult`
#' @param pcs The PC's to show (min 2, max 3). Default is `1:2`
#' @param topn the number of top genes to enumerate that drive direction of PCs
#' @param feature_id_col the column name in `x[["row_covariates"]]` to use
#'   in the feature table.
vizualize.FacilePCAResult <- function(x, pcs = 1:2, ntop = 100,
                                      with_features = TRUE,
                                      feature_id_col = NULL, ...,
                                      xlabel = "default",
                                      ylabel = "default",
                                      zlabel = "default") {
  xx <- result(x)
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

  p <- fscatterplot(xx, pc.cols, xlabel = xlabel, ylabel = ylabel,
                    zlabel = zlabel, ...)

  p$analysis <- x

  if (with_features) {
    franks <- head(ranks(x, feature_id_col, ...), ntop)
    franks <- franks[, c("rank", pc.cols)]
    dtable <- datatable(
      franks, extensions = "Scroller", style = "bootstrap",
      class = "compact", width = "100%", rownames = FALSE,
      options=list(deferRender=TRUE, scrollY=300, scroller=TRUE))
    out <- bscols(p$plot, dtable, widths = c(8, 4))
  } else {
    out <- p$plot
  }
  out
}

#' @export
#' @noRd
report.FacilePCAResult <- function(x, ...) {
  vizualize(x, ...)
}
