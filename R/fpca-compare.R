#' @rdname fpca
#' @export
#'
#' @section Comparing PCA Results:
#' We can compare two PCA results. Currently this just means we compare the
#' loadings of the features along each PC from fpca result `x` and `y`.
#'
#' @param rerun when `rerun = TRUE` (default), the `fpca(x)` and `fpca(y)` will
#'   be rerun over the union of the features in `x` and `y`.
#'
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#' p1 <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' p2 <- efds %>%
#'   FacileData::filter_samples(indication == "BLCA") %>%
#'   fpca()
#' pcmp <- compare(p1, p2)
compare.FacilePcaAnalysisResult <- function(x, y, run_all = TRUE, rerun = TRUE,
                                            ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  # TODO: assert_comparable(x, y, ...)
  # TODO: assert_comparable.FacileDgeAnalysisResult <- function(...)
  assert_class(x, "FacilePcaAnalysisResult")
  assert_class(y, "FacilePcaAnalysisResult")
  aname <- param(x, "assay_name")
  assert_true(aname == param(y, "assay_name"))

  fds. <- assert_class(fds(x), "FacileDataStore")
  assert_true(name(fds.) == name(fds(y)))

  clazz <- NULL
  classes <- c(
    "FacilePcaComparisonAnalysisResult",
    if (run_all) "FacilePcaAnalysisResult" else NULL,
    "FacileComparisonAnalysis",
    "FacileAnalysisResult")
  out <- list(
    result = NULL,
    params = list(x = x, y = y),
    fds = fds.)

  # down to busines ............................................................
  fids.x <- features(x)[["feature_id"]]
  fids.y <- features(y)[["feature_id"]]
  fids.both <- intersect(fids.x, fids.y)

  fids <- union(fids.x, fids.y)
  ndims <- max(param(x, "dims"), param(y, "dims"))
  ndims <- min(ndims, length(fids), nrow(samples(x)), nrow(samples(y)))

  if (run_all) {
    all.samples <- bind_rows(samples(x), samples(y))
    all.samples <- as_facile_frame(all.samples, fds.)
    all.samples <- distinct(all.samples, dataset, sample_id, .keep_all = TRUE)
    ures <- fpca(all.samples, assay_name = aname, dims = ndims,
                 features = fids)
    add.params <- setdiff(names(ures[["params"]]), names(out[["params"]]))

    out[["params"]] <- c(out[["params"]], ures[["params"]][add.params])
    out[["result"]] <- ures
    out[["dims"]] <- seq(ndims)
    out[["rotation"]] <- ures[["rotation"]]
    out[["percent_var"]] <- ures[["percent_var"]]
    out[["row_covariates"]] <- ures[["row_covariates"]]
    out[["takne"]] <- ures[["taken"]]
    out[["samples"]] <- ures[["samples"]]
    messages <- c(messages, ures[["messages"]])
    warnings <- c(warnings, ures[["warnings"]])
    errors <- c(errors, ures[["errors"]])
  }

  if (rerun == TRUE) {
    x <- fpca(samples(x), assay_name = aname, dims = ndims, features = fids)
    y <- fpca(samples(y), assay_name = aname, dims = ndims, features = fids)
  }

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, classes)
    return(out)
  })

  dims.take <- min(param(x, "dims"), param(y, "dims"))
  ranks.x <- tidy(ranks(x, dims = seq(dims.take)))
  ranks.y <- tidy(ranks(y, dims = seq(dims.take)))

  meta.cols <- c("symbol", "meta", "name")
  meta.cols <- intersect(meta.cols, colnames(ranks.x))
  xystats <- full_join(
    select(ranks.x, feature_id, {{meta.cols}}, dimension, rank, score, weight),
    select(ranks.y, feature_id, dimension, rank, score, weight),
    by = c("feature_id", "dimension")) %>%
    group_by(dimension) %>%
    mutate(
      membership = case_when(
        feature_id %in% fids.both ~ "both",
        feature_id %in% fids.x    ~ "only_x",
        feature_id %in% fids.y    ~ "only_y",
        TRUE                      ~ "wtf"),
      score.x = ifelse(is.na(score.x), min(score.x, na.rm = TRUE), score.x),
      score.y = ifelse(is.na(score.y), min(score.y, na.rm = TRUE), score.y)) %>%
    ungroup()


  out[["xystats"]] <- xystats
  out
}

#' @noRd
#' @export
tidy.FacilePcaComparisonAnalysisResult <- function(x, name = "xystats", ...) {
  assert_choice(name ,c("xystats", "result", "x", "y"))
  out <- switch(name,
                xystats = x[["xystats"]],
                result = tidy(result(x)),
                x = tidy(param(x, "x")),
                y = tidy(param(x, "y")))
  out
}

