#' Retrieves sample covariates with option to use custom covariates by user.
#'
#' This returns a mix of custom covs that might be provided by a passed (wide)
#' data.frame and the ones stored in the datastore. Covariates passed in by an
#' extra data.frame or the original sample descriptor are given higher
#' precedence.
#'
#' @export
#' @param x a sample descriptor (facile_frame)
#' @param covariates the names of the covariates to retrieve from the datstore.
#'   Defaults to `NULL`, which is all of them.
#' @param custom_covariates a wide (sample, dataset, cov1, cov2, ...) data.frame
#'   of extra covariates the user wants to add to the sampel descriptor `x`.
#' @param custom_key for [FacileData::with_sample_covariates()]
#' @return a wider version of `x` with more covariates attached to it.
with_sample_covs <- function(x, covariates = NULL, custom_covariates = NULL,
                             custom_key = Sys.getenv("USER"), ...,
                             verbose = FALSE) {
  assert_sample_subset(x)

  if (!is.null(custom_covariates)) {
    assert_sample_subset(custom_covariates)
    x.clash <- intersect(colnames(x), colnames(custom_covariates))
    x.clash <- setdiff(x.clash, c("dataset", "sample_id"))
    if (length(x.clash)) {
      if (verbose) {
        warning(length(x.clash), " covariates found in both the sample ",
                "descriptor and custom_covariates. ",
                "custom_covariates gets priority:\n",
                "   ", paste(x.clash, collapse = ","))
      }
      x <- select(x, -!!x.clash)
    }
    x <- left_join(x, custom_covariates, by = c("dataset", "sample_id"))
  }

  out <- with_sample_covariates(x, covariates, custom_key = custom_key)
  out
}
