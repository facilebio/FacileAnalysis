#' This is a helper function to pull out multiple signature from a ranking table
#' that provides multiple columns of rankings for the same features, like the
#' rankings produced by ranks(fpca()), where we have different rankings per PC.
#'
#' I'm not sure if this is generic enough to re-use, but I thought it might be,
#' so the signature.FacilePCAFeatureRankings function is at least using this
#' for now.
#'
#' @noRd
signature.MultiDimRankings <- function(x, ranking_columns, ntop = 20,
                                       collection_name = class(x)[1L],
                                       ...) {
  assert_class(x, "FacileFeatureRankings")
  res <- result(x)
  assert_string(collection_name)
  assert_character(ranking_columns)
  assert_subset(ranking_columns, colnames(res))

  all.rank.cols <- x[["ranking_columns"]]

  sigs <- lapply(ranking_columns, function(col) {
    out <- mutate(res, rank = res[[col]], collection = collection_name,
                  name = col)
    out <- select(out, rank, everything(), -!!all.rank.cols)
    out <- arrange(out, rank)
    head(out, ntop)
  })
  sigs <- bind_rows(sigs)
  select(sigs, collection, name, rank, everything())
}

# Print Rankings ===============================================================

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

# Print Signatures =============================================================

#' @noRd
#' @export
print.FacileFeatureSignature <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
format.FacileFeatureSignature <- function(x, ...) {
  out <- paste(
    "===========================================================\n",
    class(x)[1L], "\n",
    "-----------------------------------------------------------\n",
    "Number of features: ", nrow(result(x)), "\n",
    "===========================================================\n",
    sep = "")
  out
}
