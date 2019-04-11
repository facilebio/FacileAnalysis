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
#' @importFrom multiGSEA eigenWeightedMean
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
fpca <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                 col_covariates = NULL, ...) {
  UseMethod("fpca", x)
}

#' @noRd
#' @export
fpca.FacileDataStore <- function(x, pcs = 1:10, ntop = 500,
                                 row_covariates = NULL, col_covariates = NULL,
                                 assay_name = default_assay(x),
                                 custom_key = Sys.getenv("USER"), ...) {
  fpca(samples(x), pcs, ntop, row_covariates, col_covariates, assay_name,
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
fpca.facile_frame <- function(x, pcs = 1:10, ntop = 500,
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
    select(dataset, sample_id, everything())
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
fpca.DGEList <- function(x, pcs = 1:10, ntop = 500, row_covariates = x$genes,
                         col_covariates = x$samples,
                         prior.count = 3, assay_name = "counts", ...) {
  if (assay_name == "counts") {
    m <- edgeR::cpm(x, prior.count = prior.count, log = TRUE)
  } else {
    m <- x[[assay_name]]
    assert_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))
  }
  out <- fpca(m, pcs, ntop, row_covariates, col_covariates, ...)
}

#' @export
#' @rdname fpca
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, pcs = 1:10, ntop = 500, row_covariates = NULL,
                        col_covariates = NULL, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

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

  factor_contrib <- ewm$factor.contrib %>%
    rename(feature_id = "featureId") %>%
    select(feature_id, !!pcs.take) %>%
    as.tbl()

  result <- list(
    result = dat,
    pcs = pcs.take,
    factor_contrib = factor_contrib,
    percent_var = percentVar,
    pc_cor = pc_cor,
    row_covariates = row_covariates,
    taken = take,
    # Standard FacileAnalysisResult things
    # fds = .fds,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(result) <- c("FacilePCAResult", "FacileReducedDimResult",
                     "FacileAnalysisResult")
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
ranks.FacilePCAResult <- function(x, type = c("rankings", "ranked"),
                                  report_feature_as = NULL,  ...) {
  type <- match.arg(type)
  fcontrib <- x[["factor_contrib"]]
  rdata <- x[["row_covariates"]]

  # pc.cols <- colnames(fcontrib)[grepl("PC\\d+", colnames(fcontrib))]
  pc.cols <- x[["pcs"]]
  pc.ranks <- fcontrib[, c("feature_id", pc.cols)]
  for (pc in pc.cols) {
    pc.ranks[[pc]] <- rank(-pc.ranks[[pc]], ties.method = "random")
  }

  if (!"feature_id" %in% colnames(rdata)) {
    rdata[["feature_id"]] <- rownames(fcontrib)
  }
  if (!"feature_type" %in% colnames(rdata)) {
    # TODO: extract the feature_type of these rankings from the
    # "facile analysis chain"
    feature_type <- guess_feature_type(rdata[["feature_id"]], summarize = TRUE)
    rdata[["feature_type"]] <- feature_type[["feature_type"]]
  }

  # rankings <- as.tbl(left_join(pc.ranks, rdata, by = "feature_id"))
  rankings <- pc.ranks %>%
    left_join(rdata, by = "feature_id") %>%
    select(feature_type, feature_id, everything())

  if (type == "rankings") {
    result. <- rankings
  } else if (type == "ranked") {
    meta.cols <- setdiff(colnames(rankings), pc.cols)
    if (is.null(report_feature_as)) {
      opts <- c("name", "symbol", "feature_id", "featureId",
                "gene_id", "geneId")
      report_feature_as <- intersect(opts, meta.cols)[1L]
      if (is.na(report_feature_as)) {
        report_feature_as <- meta.cols[1L]
      }
    }
    assert_choice(report_feature_as, meta.cols)
    result. <- tibble(rank = seq(nrow(rankings)))
    for (pc in pc.cols) {
      ordr <- order(rankings[[pc]])
      result.[[pc]] <- rankings[[report_feature_as]][ordr]
    }
  }

  out <- list(
    result = as.tbl(result.),
    feature_weight = gather(fcontrib, PC, weight, -feature_id),
    # TODO: add faciledatastore to ranks.fpca output?
    ranking_columns = pc.cols,
    ranking_order = "ascending",
    percent_var = x[["percent_var"]])

  clazz <- paste0("Facile%sFeature", tools::toTitleCase(type))
  classes <- sprintf(clazz, c("PCA", ""))
  class(out) <- c(classes, "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
signature.FacilePCAFeatureRankings <- function(x, pcs = NULL, ntop = 20,
                                               collection_name = class(x)[1L],
                                               ranking_columns = x[["ranking_columns"]],
                                               ...) {
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

  sigs <- signature.MultiDimRankings(x, ranking_columns, ntop = ntop,
                                     collection_name = collection_name, ...)
  # Add PCA specific metadata to features
  pvar <- x[["percent_var"]]
  sigs <- left_join(sigs, x[["feature_weight"]],
                    by = c("name" = "PC", "feature_id"))
  sigs[["percent_var"]] <- pvar[sigs[["name"]]]
  front.load <- c("collection", "name", "feature_id", "symbol", "rank",
                  "weight", "percent_var")
  front.load <- intersect(front.load, colnames(sigs))
  sigs <- select(sigs, !!front.load, everything())
  out <- list(
    result = sigs,
    params = list(pcs = pcs, ntop = ntop))
  class(out) <- c("FacilePCAFeatureSignature",
                  "FacileFeatureSignature",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
signature.FacilePCAResult <- function(x, pcs = NULL, ntop = 20, ...) {
  signature(ranks(x, "rankings"), pcs = pcs, ntop = ntop, ...)
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

