# Most analyses can induce well defined rankings over the feature space that
# is analyzed. These are methods that work over such rankings, ie.
#
# `ranks(x, type = "rankings")`

#' Convert a feature descriptor into a GeneSetDb
#'
#' This is mainly for feature rankings, but could be something else?
#'
#' @export
#' @examples
#' pca.gdb <- FacileData::exampleFacileDataSet() %>%
#'   filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   ranks() %>%
#'   as.GeneSetDb()
as.GeneSetDb <- function(x, ...) {
  UseMethod("as.GeneSetDb", x)
}

#' @export
#' @importFrom multiGSEA GeneSetDb
#' @noRd
#' @method as.GeneSetDb FacileFeatureRankings
as.GeneSetDb.FacileFeatureRankings <- function(x, topn = 10,
                                               collection_name = class(x)[1L],
                                               ranking_columns = x[["ranking_columns"]],
                                               ...) {
  res <- result(x)
  assert_string(collection_name)
  assert_character(ranking_columns)
  assert_subset(ranking_columns, colnames(res))

  all.rank.cols <- x[["ranking_columns"]]

  sigs <- lapply(ranking_columns, function(col) {
    out <- mutate(res, rank = res[[col]], collection = collection_name,
                  name = col)
    out <- select(out, rank, everything(), -!!all.rank.cols)
    out <- rename(out, featureId = "feature_id")
    out <- arrange(out, rank)
    head(out, topn)
  })
  sigs <- bind_rows(sigs)
  GeneSetDb(sigs)
}

#' @noRd
#' @export
#' @method as.GeneSetDb FacilePCAFeatureRankings
as.GeneSetDb.FacilePCAFeatureRankings <- function(x, pcs = NULL, ...) {
  res <- result(x)
  if (is.null(pcs)) {
    ranking_columns <-  x[["ranking_columns"]]
  } else if (test_int(pcs)) {
    ranking_columns <- paste0("PC", 1:pcs)
  } else if (test_integerish(pcs)) {
    ranking_columns <- paste0("PC", pcs)
  }
  assert_character(ranking_columns)
  assert_subset(ranking_columns, colnames(res))
  gdb <- NextMethod(x, ranking_columns = ranking_columns, ...)

  # Let's add the percent_var explained for each PC to the metadata for the
  # genesets
  pvar <- x[["percent_var"]]
  gdb@table[["percent_var"]] <- pvar[gdb@table[["name"]]]
  gdb
}

#' @noRd
#' @export
print.FacileFeatureRankings <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacileFeatureRankings <- function(x, ...) {
  out <- paste(
    "===========================================================\n",
    class(x)[1L], "\n",
    "-----------------------------------------------------------\n",
    "Number of features: ", nrow(result(x)), "\n",
    "===========================================================\n",
    sep = "")
  out
}
