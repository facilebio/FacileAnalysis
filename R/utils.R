# A place for random utility functions

#' Extract a feature_id character vector from a feature descriptor, or NULL
#'
#' All of the functions that take in a `features` parameter can call this
#' function to do the id extraction ... I wrote this code a lot and all
#' over the place.
#'
#' @noRd
#' @param x NULL, character string, or tibble with feature_id column
#' @return a character vector of feature ids
extract_feature_id <- function(x, as_tibble = FALSE, ...) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- x[["feature_id"]]
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) {
    stop("Could not extract feature_id from feature descriptor", call. = FALSE)
  }
  if (as_tibble) {
    x <- tibble(feature_id = x)
  }
  x
}

#' The primary key for a feature is ensembl id, uniprotid, etc, but somtimes you want to look up by gene names
#' 
#' @export
#' @param x a character list of feature names/symbols
#' @param reference an object you are looking up in, like a FacileAnalysisResult, facile_frame, fds, ...
#' @export
#' @examples
#' feature_lookup(c("MCM2", "ENSG00000088930"), FacileData::an_fds())
feature_lookup <- function(x, reference, assay_name = NULL, feature_name = "name", ...) {
  if (is.character(x)) {
    x <- dplyr::tibble(feature_id = x)
  }
  stopifnot(
    "data.frame required" = is.data.frame(x),
    "feature_id column required" = is.character(x[["feature_id"]])
  )

  # dispatch on reference (sigh ...)
  fids.all <- NULL
  fds. <- NULL
  if (is(reference, "FacileAnalysisResult")) {
    fds. <- fds(reference)
  } else if (is(reference, "facile_frame")) {
    fds. <- fds(reference)
  } else if (is(reference, "FacileDataStore")) {
    fds. <- reference
  }
  
  if (is.null(fds.)) {
    # user sent in a data.frame of features to search against
    fids.all <- reference
  } else {
    fids.all <- FacileData::features(fds(reference), assay_name = assay_name)  
  }
  
  checkmate::assert_data_frame(fids.all)
  checkmate::assert_subset(c(feature_name, "feature_id"), colnames(fids.all))
  fids.all <- fids.all |>
    dplyr::mutate(
      match_name = tolower(.data[[feature_name]]),
      match_id = tolower(feature_id),
      feature_row_index = 1:nrow(fids.all)
    )
  fids.all$feature_name <- fids.all[[feature_name]]
  
  x <- x |> 
    dplyr::mutate(flower. = tolower(feature_id)) |> 
    dplyr::rename(fid. = feature_id)
  
  fid.out <- x |> 
    dplyr::inner_join(
      fids.all,
      by = c("flower." = "match_id"),
      suffix = c("", ".dropme")) |> 
    dplyr::select(-dplyr::ends_with(".dropme"))

  name.out <- x |> 
    dplyr::inner_join(
      fids.all,
      by = c("flower." = "match_name"),
      suffix = c("", ".dropme")) |> 
    dplyr::select(-dplyr::ends_with(".dropme"))
  # 
  # ens.id <- grepl("^ENS[A-Z]*G\\d+$", x$feature_id, ignore.case = TRUE)
  # not.ens <- x |>
  #   dplyr::filter(!ens.id) |>
  #   dplyr::rename(match_name = feature_id) |>
  #   dplyr::inner_join(
  #     fids.all,
  #     by = "match_name",
  #     suffix = c("", ".dropme")) |> 
  #   dplyr::select(-dplyr::ends_with(".dropme"))
  # 
  # ens <- x |>
  #   dplyr::filter(ens.id) |>
  #   dplyr::inner_join(
  #     fids.all,
  #     by = c("feature_id" = "match_id"),
  #     suffix = c("", ".dropme")
  #   ) |>
  #   dplyr::mutate(feature_id = feature_id.dropme) |> 
  #   dplyr::select(-dplyr::ends_with(".dropme"))
  # 
  # out <- not.ens |>
  #   dplyr::bind_rows(ens) |>
  #   dplyr::arrange(name) |>
  #   dplyr::select(!dplyr::starts_with("match_"))
  
  out <- fid.out |> 
    bind_rows(name.out) |> 
    dplyr::select(!dplyr::starts_with("match_")) |> 
    dplyr::rename(feature_query = fid.) |> 
    dplyr::relocate(feature_query, .after = dplyr::last_col())
  
  not.found <- !x$flower. %in% c(tolower(out$feature_id), tolower(out[[feature_name]]))
  if (any(not.found)) {
    warning(
      "These features were not found: ",
      paste(x$fid.[not.found], collapse = ",")
    )
  }
  
  out <- dplyr::select(out, -flower.)
  attr(out, "missing") <- dplyr::select(x[not.found,,drop=FALSE], query = fid.)
  out
}


#' Convenience wrapper to require specified packages
#'
#' @noRd
#' @param pkg A character vector of packages to require
#' @param quietly defaults to true
#' @param ... passed into [requireNamespace()]
reqpkg <- function(pkg, quietly = TRUE, ...) {
  assert_character(pkg)
  for (p in pkg) {
    if (!requireNamespace(p, ..., quietly = quietly)) {
      stop("'", p, "' package required, please install it.", call. = FALSE)
    }
  }
}
