.reduce_dim_methods <- function(x) {
  c("pca", "tsne", "nmf")
}

#' General dimensionality reduction
#'
#' @export
#' @rdname freduce_dims
#' @return a `FacileReducedDimsResult`
freduce_dims <- function(x, method = .reduce_dim_methods(), ...) {
  UseMethod("freduce_dims", x)
}
