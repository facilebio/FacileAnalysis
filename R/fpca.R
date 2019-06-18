#' Runs a facile PCA
#'
#' The code here is largely inspired by DESeq2's plotPCA.
#'
#' You should look at factominer:
#' * http://factominer.free.fr/factomethods/index.html
#' * http://factominer.free.fr/graphs/factoshiny.html
#'
#' @section Teaching and Tutorials:
#'
#' This looks like a useful tutorial to use when explaining the utility of
#' PCA analysis:
#' http://alexhwilliams.info/itsneuronalblog/2016/03/27/pca/
#'
#' @export
#' @rdname fpca
#'
#' @param x a data container
#' @return an fpca result
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#'
#' # A subset of samples ------------------------------------------------------
#' pca.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' if (interactive()) {
#'   report(pca.crc, color_aes = "sample_type")
#' }
#'
#' pca.gdb <- pca.crc %>%
#'   signature(pcs = 1:3) %>%
#'   result() %>%
#'   multiGSEA::GeneSetDb()
#'
#' # All samples --------------------------------------------------------------
#' pca.all <- fpca(efds)
#' if (interactive()) {
#'   viz(pca.all, color_aes = "indication", shape_aes = "sample_type")
#'   report(pca.all, color_aes = "indication", shape_aes = "sample_type")
#' }
#'
#'
#' # This works on "normal" DGELists, too. -----------------------------------
#' pca.dgelist <- efds %>%
#'   filter_samples(indication == "CRC") %>%
#'   as.DGEList() %>%
#'   fpca()
#' if (interactive()) {
#'   report(pca.dgelist, color_aes = "sample_type")
#' }
fpca <- function(x, pcs = 5, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("fpca", x)
}

#' @noRd
#' @export
fpca.FacileDataStore <- function(x, pcs = 5, ntop = 500,
                                 row_covariates = NULL, col_covariates = NULL,
                                 assay_name = default_assay(x),
                                 custom_key = Sys.getenv("USER"), ...,
                                 samples = NULL) {
  if (is.null(samples)) samples <- samples(x)
  samples <- collect(samples, n = Inf)

  fpca(samples, pcs, ntop, row_covariates, col_covariates, assay_name,
       custom_key, ...)
}

#' @section FacileDataStore (facile_frame):
#' We enable the user to supply extra sample covariates that are not found
#' in the FacileDataStore associated with these samples `x` by adding them as
#' extra columns to `x`.
#'
#' If manually provioded col_covariates have the same name as internal sample
#' covariates, then the manually provided ones will supersede the internals.
#'
#' @rdname fpca
#' @export
fpca.facile_frame <- function(x, pcs = 5, ntop = 500,
                              row_covariates = NULL, col_covariates = NULL,
                              assay_name = NULL,
                              custom_key = Sys.getenv("USER"), ...) {
  .fds <- assert_class(fds(x), "FacileDataStore")
  assert_sample_subset(x)
  x <- collect(x, n = Inf)

  if (!is.null(row_covariates)) {
    warning("Custom row_covariates not yet supported for facile_frame ",
            "(it's not hard, I'm just lazy right now", immediate. = TRUE)
  }

  # TODO: Now that the result() of this thing is a facile_frame, I don't think
  # we need to add these covariates
  col.covariates <- with_sample_covs(x, custom_covariates = col_covariates,
                                     custom_key = custom_key)

  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
  }

  # (lazily) turning this into a DGEList to leverage the already-implented
  # feature and sample anntation alignment written up in there.
  y <- as.DGEList(x, covariates = col.covariates, assay_name = assay_name)

  out <- fpca(y, pcs, ntop, ...)
  out[["result"]] <- out[["result"]] %>%
    as.tbl() %>%
    select(dataset, sample_id, everything()) %>%
    as_facile_frame(.fds)
  out[["feature_stats"]] <- out[["feature_stats"]] %>%
    as.tbl() %>%
    as_facile_frame(.fds)

  out[["params"]] <- list(
    pcs = pcs,
    ntop = ntop,
    assay_name = assay_name)
  ## add facile stuff
  out[["fds"]] <- .fds
  out
}

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
fpca.DGEList <- function(x, pcs = 5, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = 3, assay_name = "counts", ...) {
  if (assay_name == "counts") {
    m <- edgeR::cpm(x, prior.count = prior.count, log = TRUE)
  } else {
    m <- x[[assay_name]]
    assert_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))
  }
  out <- fpca(m, pcs, ntop, row_covariates, col_covariates, ...)
  out
}

#' @noRd
#' @export
fpca.EList <- function(x, pcs = 5, ntop = 500, row_covariates = x$genes,
                       col_covariates = x$targets,
                       prior.count = 3, assay_name = "counts", ...) {
  out <- fpca(x$E, pcs, ntop, row_covariates, col_covariates, ...)
  out
}

#' @export
#' @rdname fpca
#' @importFrom irlba prcomp_irlba
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, pcs = 5, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL, use_irlba = pcs < 7,
                        center = TRUE, scale. = FALSE, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  assert_int(pcs, lower = 3L, upper = min(nrow(x), ncol(x)))
  assert_flag(use_irlba)

  if (is.null(rownames(x))) rownames(x) <- as.character(seq(nrow(x)))
  if (is.null(row_covariates)) {
    row_covariates <- data.frame(feature_id = rownames(x),
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

  rv <- matrixStats::rowVars(x)
  take <- head(order(rv, decreasing = TRUE), ntop)

  xx <- x[take,,drop = FALSE]
  row_covariates <- row_covariates[take,,drop = FALSE]

  if (use_irlba) {
    pca <- prcomp_irlba(t(xx), n = pcs, center = center, scale. = scale.)
    rownames(pca$rotation) <- rownames(xx)
    rownames(pca$x) <- colnames(xx)
  } else {
    pca <- prcomp(t(xx), center = center, scale. = scale.)
    pca$sdev <- head(pca$sdev, pcs)
    pca$rotation <- pca$rotation[, 1:pcs, drop = FALSE]
    pca$x <- pca$x[, 1:pcs, drop = FALSE]
  }

  pca$sdev <- setNames(pca$sdev, paste0("PC", seq(pca$sdev)))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  dat <- as.data.frame(pca$x)

  if (is(col_covariates, "data.frame")) {
    dat <- cbind(dat, col_covariates[rownames(dat),,drop = FALSE])
  }

  result <- list(
    result = dat,
    pcs = pcs,
    rotation = pca$rotation,
    percent_var = percentVar,
    row_covariates = row_covariates,
    taken = take,
    # Standard FacileAnalysisResult things
    # fds = .fds,
    messages = messages,
    warnings = warnings,
    errors = errors)
  class(result) <- c("FacilePCAResult", "FacileReducedDimResult",
                     "FacileAnalysisResult")

  result[["feature_stats"]] <- .fpca.feature_statistics(result)
  result
}

# Ranks and Signatures =========================================================

#' Reports the contribution of each gene to the principal components
#'
#' NOTE: type = "ranked" shouldn't be here. It is used by the viz
#' display the genes loaded on the PCs in a wide format for visual
#' function. We need to think of a universal function name that will do
#' this in the "general" sense.
#'
#' @export
#' @noRd
#' @param x A `FacilePCAResult`
#' @param type `"rankings"` (default) or `"ranked"` -- note that the idea
#' of the output for `"ranked"` should be generalized. See `NOTE` in the
#' function details section.
#' @param report_feature_as The column used to display in the returned
#'   table for each feature when `type == "rankded"`
#' @return A FacilePCAFeature(Rankings|Ranked) object
ranks.FacilePCAResult <- function(x, type = c("features", "samples"),
                                  signed = TRUE, ...) {
  type <- match.arg(type)
  if (type == "samples") stop("What does sample-ranking even mean?")

  if (type == "features") {
    fstats <- x[["feature_stats"]]
    rcol <- if (signed) "rank_rotation" else "rank_weight"
    ranks. <- select(fstats, feature_id, feature_type,
                     dimension = PC, score = rotation,
                     weight, rank = !!rcol)
    ranks. <- arrange(ranks., dimension, rank)
    clazz <- "FacilePCAFeatureRankings"
  }

  pvar <- x[["percent_var"]]
  ranks.[["percent_var"]] <- pvar[ranks.[["dimension"]]]

  out <- list(
    result = ranks.,
    params = list(type = type, signed = signed),
    fds = fds(x))

  class(out) <- c(clazz,
                  "FacileMultiDimensionalRankings",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
signature.FacilePCAFeatureRankings <- function(x, pcs = NULL, ntop = 20,
                                               collection_name = class(x)[1L],
                                               ranking_columns = x[["ranking_columns"]],
                                               ...) {
  if (isTRUE(x$params$signed)) {
    sig.up <- result(x) %>%
      group_by(dimension) %>%
      arrange(rank) %>%
      slice(1:min(ntop, n())) %>%
      mutate(name = paste(dimension, "pos")) %>%
      ungroup()
    sig.down <- result(x) %>%
      group_by(dimension) %>%
      arrange(desc(rank)) %>%
      slice(1:min(ntop, n())) %>%
      mutate(name = paste(dimension, "neg")) %>%
      ungroup()
    sig <- sig.up %>%
      bind_rows(sig.down) %>%
      arrange(dimension, rank)
  } else {
    sig <- result(x) %>%
      group_by(dimension) %>%
      arrange(rank) %>%
      slice(1:min(ntop, n())) %>%
      mutate(name = paste(dimension, "unsigned")) %>%
      ungroup()
  }

  sig <- mutate(sig, collection = collection_name)
  sig <- select(sig, collection, name, feature_id, everything())
  sig <- as_facile_frame(sig, fds(x))

  out <- list(
    result = sig,
    params = list(pcs = pcs, ntop = ntop))

  # FacilePCAFeatureSignatureSigned
  clazz <- "Facile%sFeatureSignature%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c("PCA", "PCA",  "", ""), c(s, "", s, ""))
  class(out) <- c(classes,
                  "FacileFeatureRanks",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
signature.FacilePCAResult <- function(x, type = "features", signed = TRUE,
                                      pcs = NULL, ntop = 20, ...) {
  signature(ranks(x, type = type, signed = signed), pcs = pcs, ntop = ntop, ...)
}

# Utility Functions ============================================================

#' Helper function that creates a long/tidy table of statistics over the
#' features of a PCA (loadings, weights, ranks)
#' @noRd
.fpca.feature_statistics <- function(x, ...) {
  assert_class(x, "FacilePCAResult")

  pvar <- assert_numeric(x$percent_var, names = "unique")
  pc.names <- names(pvar)
  rotation <- expect_matrix(x[["rotation"]], mode = "numeric",
                            ncols = length(pvar), col.names = "unique")
  assert_subset(pc.names, colnames(rotation))
  rotation <- bind_cols(tibble(feature_id = rownames(rotation)),
                        as.data.frame(rotation))

  rdata <- x[["row_covariates"]]
  if (!"feature_id" %in% colnames(rdata)) {
    rdata[["feature_id"]] <- rownames(rdata)
  }
  assert_character(rdata[["feature_id"]], unique = TRUE)
  assert_true(all(rdata[["feature_id"]] == rotation[["feature_id"]]))


  meta.cols <- c("feature_id", "feature_type")
  meta.cols <- intersect(meta.cols, colnames(rdata))
  assert_character(meta.cols, min.len = 1L)
  meta <- rdata[, meta.cols, drop = FALSE]
  assert_character(meta[["feature_id"]])

  if (!is.character(meta[["feature_type"]])) {
    ftype <- guess_feature_type(meta[["feature_id"]],
                                with_organism = FALSE,
                                summarize = TRUE)
    meta[["feature_type"]] <- ftype[["feature_type"]]
  }

  rlong <- rotation %>%
    gather("PC", "rotation", -feature_id) %>%
    mutate(weight = abs(rotation)) %>%
    group_by(PC) %>%
    mutate(rank_rotation = rank(-rotation, ties.method = "random"),
           rank_weight = rank(-weight, ties.method = "random")) %>%
    ungroup()

  stats <- inner_join(meta, rlong, by = "feature_id")
  as.tbl(stats)
}

# Printing =====================================================================

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
    "Number of features used: ", n.features, "\n",
    "Number of PCs: ", length(pcv), "\n",
    "Variance explained:\n  ", pcvu, "\n",
    "===========================================================\n",
    sep = "")
  out
}

