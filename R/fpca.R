#' Runs a principal components analysis, the facile way.
#'
#' Performs a principal components analysis over a specified assay from the
#' (subset of) samples in a FacileDataStore.
#'
#' The `FacilePcaAnalysisResult` produced here can be used in "the usual" ways,
#' ie. can be `viz`-ualized. `shine()` is 1/4th-implemented, and `report()`
#' has not been worked on yet.
#'
#' Importantly / interestingly, you can shoot this result into [ffsea()] to
#' perform gene set enrichment analysis over a specified dimension to identify
#' functional categories loaded onto differend PCs.
#'
#' @section Development Notes:
#' Follow progress on implementation of `shine()` and `report()` below:
#'
#' 1. [Implement `shine()`](https://github.com/facileverse/FacileAnalysis/issues/11)
#' 2. [Implement `report()`](https://github.com/facileverse/FacileAnalysis/issues/12)
#'
#' Note that there are methods defined for other assay containers, like an
#' `edgeR::DGEList`, `limma::EList`, and `SummarizedExperiment`. If these are
#' called directly, their downstream use within the facile ecosystem isn't
#' yet fully supported. Development of the
#' [FacileBioc package](https://github.com/facileverse/FacileBioc)
#' will address this.
#'
#' @section Random Things to elaborate on:
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
#' @param x a facile data container (FacileDataSet), or a `facile_frame`
#'   (refer to the FacileDataStore (facile_frame) section.
#' @param assay_name the name of the assay to extract data from to perform the
#'   PCA. If not specified, default assays are taken for each type of assay
#'   container (ie. `default_assay(facile container)`, `"counts"` for a
#'   `DGEList`, `assayNames(SummarizedExperiment)[1L]`, etc.)
#' @param dims the number of PC's to calculate (minimum is 3).
#' @param ntop the number of features (genes) to include in the PCA. Genes are
#'   ranked by decreasing variance across the samples in `x`.
#' @param row_covariates,col_covariates data.frames that provie meta information
#'   for the features (rows) and samples (columns). The default is to get
#'   these values from "the obvious places" given `x` (`$genes` and `$samples`
#'   for a DGEList, or the sample and feature-level covariate database tables
#'   from a FacileDataSet, for example).
#' @param batch,main specify the covariates to use for batch effect removal.
#'   Refer to the [FacileData::remove_batch_effect()] help for more information.
#' @return an fpca result
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#'
#' # A subset of samples ------------------------------------------------------
#' pca.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' if (interactive()) {
#'   # report(pca.crc, color_aes = "sample_type")
#'   shine(pca.crc)
#'   viz(pca.crc, color_aes = "sex")
#' }
#'
#' # Same PCA as above, but regress "sex" out of samples first.
#' pca.crcs <- samples(pca.crc) %>%
#'   fpca(batch = "sex")
#' if (interactive()) {
#'   viz(pca.crcs, color_aes = "sex")
#' }
#'
#' pca.gdb <- pca.crc %>%
#'   signature(dims = 1:3) %>%
#'   result() %>%
#'   multiGSEA::GeneSetDb()
#'
#' # All samples --------------------------------------------------------------
#' pca.all <- fpca(efds)
#' if (interactive()) {
#'   viz(pca.all, color_aes = "indication", shape_aes = "sample_type")
#'   # report(pca.all, color_aes = "indication", shape_aes = "sample_type")
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
fpca <- function(x, assay_name = NULL, dims = 5, filter = "default",
                 ntop = 1000, row_covariates = NULL, col_covariates = NULL,
                 batch = NULL, main = NULL, ...) {
  UseMethod("fpca", x)
}

#' @noRd
#' @export
fpca.FacileDataStore <- function(x, assay_name = NULL, dims = 5,
                                 filter = "default", ntop = 1000,
                                 row_covariates = NULL,
                                 col_covariates = NULL, batch = NULL,
                                 main = NULL, custom_key = Sys.getenv("USER"),
                                 ..., samples = NULL) {
  assert_int(dims, lower = 3L) # TODO: max = min(dim(x, assay_name = ??))
  if (is.null(samples)) samples <- samples(x)
  samples <- collect(samples, n = Inf)

  fpca(samples, assay_name, dims, filter, ntop, row_covariates, col_covariates,
       batch, main, custom_key, ...)
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
fpca.facile_frame <- function(x, assay_name = NULL,
                              dims = min(5, nrow(x) - 1L),
                              filter = "default", ntop = 1000,
                              row_covariates = NULL,
                              col_covariates = NULL, batch = NULL, main = NULL,
                              custom_key = Sys.getenv("USER"), ...) {
  assert_int(dims, lower = 3L) # TODO: max = min(dim(x, assay_name = ??))
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
  # (lazily) turning this into a DGEList to leverage the already-implemented
  # feature and sample anntation alignment written up in there. The desired
  # assay data (from assay_matrix) will be stored in the DGEList's $count
  # element
  y <- as.DGEList(x, covariates = col.covariates, assay_name = assay_name)

  out <- fpca(y, dims = dims, filter = filter, ntop = ntop, batch = batch,
              main = main, ...)
  out[["params"]][["assay_name"]] <- assay_name

  out[["result"]] <- out[["result"]] %>%
    as.tbl() %>%
    select(dataset, sample_id, everything()) %>%
    as_facile_frame(.fds)
  out[["feature_stats"]] <- out[["feature_stats"]] %>%
    as.tbl() %>%
    as_facile_frame(.fds)

  out[["samples"]] <- x

  # add facile stuff
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

#' @export
#' @rdname fpca
#' @importFrom irlba prcomp_irlba
#' @importFrom matrixStats rowVars
fpca.matrix <- function(x, dims = min(5, ncol(x) - 1L),
                        filter = "default", ntop = 1000,
                        row_covariates = NULL, col_covariates = NULL,
                        batch = NULL, main = NULL, use_irlba = dims < 7,
                        center = TRUE, scale. = FALSE, ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  if (min(dim(x)) < 2L) stop("Can't run PCA on a one-dimensional matrix")
  # When using irlba, n has to be strictly less than min(dim(xx))
  assert_int(dims, lower = 1L, upper = min(dim(x)) - 1L)
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

  if (!is.null(batch)) {
    x <- remove_batch_effect(x, col_covariates, batch = batch, ...)
  }

  if (test_string(filter) && filter == "default") {
    rv <- matrixStats::rowVars(x)
    take <- head(order(rv, decreasing = TRUE), ntop)
  } else {
    if (test_multi_class(filter, c("data.frame", "tbl"))) {
      filter <- filter[["feature_id"]]
    }
    assert_character(filter, min.len = 2)
    take <- match(filter, rownames(x))
    if (any(is.na(take))) {
      stop("The pca filtering strategy only allows you to specify rownames ",
           "(features) to use for PCA, and the ones you specified do not ",
           "exactly match rownames(x), see:\n  ",
           "https://github.com/facileverse/FacileAnalysis/issues/20")
    }
    ntop <- length(take)
  }

  xx <- x[take,,drop = FALSE]
  row_covariates <- row_covariates[take,,drop = FALSE]

  if (use_irlba) {
    pca <- prcomp_irlba(t(xx), n = dims, center = center, scale. = scale.)
    rownames(pca$rotation) <- rownames(xx)
    rownames(pca$x) <- colnames(xx)
  } else {
    pca <- prcomp(t(xx), center = center, scale. = scale.)
    pca$sdev <- head(pca$sdev, dims)
    pca$rotation <- pca$rotation[, 1:dims, drop = FALSE]
    pca$x <- pca$x[, 1:dims, drop = FALSE]
  }

  pca$sdev <- setNames(pca$sdev, paste0("PC", seq(pca$sdev)))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  dat <- as.data.frame(pca$x)

  if (is(col_covariates, "data.frame")) {
    dat <- cbind(dat, col_covariates[rownames(dat),,drop = FALSE])
  }
  if (!"sample_id" %in% colnames(dat)) {
    dat[["sample_id"]] <- rownames(dat)
  }

  samples. <- tibble(dataset = "dataset", sample_id = colnames(x))

  result <- list(
    result = dat,
    dims = seq(dims),
    rotation = pca$rotation,
    percent_var = percentVar,
    row_covariates = row_covariates,
    taken = take,
    samples = samples.,
    # Standard FacileAnalysisResult things
    # fds = .fds,
    params = list(dims = dims, filter = filter, ntop = ntop,
                  use_irlba = use_irlba,
                  batch = batch, main = main,
                  center = center, scale. = scale.),
    messages = messages,
    warnings = warnings,
    errors = errors)
  class(result) <- c("FacilePcaAnalysisResult",
                     "FacileReducedDimAnalysisResult",
                     "FacileAnalysisResult")

  result[["feature_stats"]] <- .fpca.feature_statistics(result)
  result
}

# Methods and Accessors ========================================================

#' @noRd
#' @export
initialized.FacilePcaAnalysisResult <- function(x, ...) {
  stat.table <- tidy(x)
  is.data.frame(stat.table) &&
    nrow(stat.table) == nrow(samples(x))
}

#' @noRd
#' @export
features.FacilePcaAnalysisResult <- function(x, ...) {
  out <- assert_class(x[["feature_stats"]], "data.frame")
  out
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
#' @param x A `FacilePcaAnalysisResult`
#' @param type `"rankings"` (default) or `"ranked"` -- note that the idea
#' of the output for `"ranked"` should be generalized. See `NOTE` in the
#' function details section.
#' @param report_feature_as The column used to display in the returned
#'   table for each feature when `type == "rankded"`
#' @return A FacilePCAFeature(Rankings|Ranked) object
ranks.FacilePcaAnalysisResult <- function(x, type = c("features", "samples"),
                                          signed = TRUE, dims = x[["dims"]][1L],
                                          ...) {
  type <- match.arg(type)
  if (type == "samples") stop("What does sample-ranking even mean?")
  if (!is.null(dims)) {
    if (is.character(dims)) { # We can accept "PC1" or 1
      dims <- as.integer(sub("^PC", "", dims))
    }
    dims <- assert_integerish(dims, lower = 1, any.missing = FALSE)
    dims <- unique(dims)
  }

  if (type == "features") {
    fstats <- x[["feature_stats"]]
    rcol <- if (signed) "rank_rotation" else "rank_weight"
    ranks. <- select(fstats, feature_id, feature_type,
                     dimension = PC, score = rotation,
                     weight, rank = {{rcol}})
    # FacilePcaFeatureRanksSigned
    clazz <- "FacilePcaFeatureRanks%s"
    s <- if (signed) "Signed" else "Unsigned"
    classes <- sprintf(clazz, c(s, ""))
    classes <- c(classes,
                 sub("Pca", "MultiDimensional", classes),
                 "FacileFeatureRanks")
  }

  if (!is.null(dims)) {
    ranks. <- filter(ranks., dimension %in% sprintf("PC%d", dims))
    if (nrow(ranks.) == 0L) {
      stop("All PC dimensions have been filtered out.")
    }
  }

  # Order by PC and rank
  ranks. <- ranks. %>%
    mutate(PC. = as.integer(sub("PC", "", ranks.$dimension))) %>%
    arrange(PC., rank) %>%
    mutate(PC. = NULL)

  # Add metadata to ranks, if there.
  rcovs <- x[["row_covariates"]]
  add.meta <- c("feature_id", setdiff(colnames(rcovs), colnames(ranks.)))
  if (length(add.meta)) {
    ranks. <- left_join(ranks., rcovs[, add.meta, drop = FALSE],
                        by = "feature_id")
  }
  pvar <- x[["percent_var"]]
  ranks.[["percent_var_dim"]] <- pvar[ranks.[["dimension"]]]

  # Put informative columns up front
  upfront <- c("feature_id", "name", "symbol", "dimension", "rank", "score",
               "weight")
  upfront <- intersect(upfront, colnames(ranks.))
  ranks. <- select(ranks., {{upfront}}, everything())

  out <- list(
    result = ranks.,
    signed = signed,
    params = list(type = type, signed = signed),
    ranking_columns = rcol,
    ranking_order = "descending",
    fds = fds(x))

  class(out) <- classes
  out
}

#' @noRd
#' @export
signature.FacilePcaFeatureRanks <- function(x, dims = NULL, ntop = 20,
                                            collection_name = class(x)[1L],
                                            ranking_columns = x[["ranking_columns"]],
                                            ...) {
  res. <- tidy(x)

  if (!is.null(dims)) {
    dims <- assert_integerish(dims, lower = 1)
    dims <- unique(dims)
    res. <- filter(res., dimension %in% sprintf("PC%d", dims))
    if (nrow(ranks.) == 0L) {
      stop("All PC dimensions have been filtered out.")
    }
  }

  if (isTRUE(x$params$signed)) {
    sig.up <- res. %>%
      group_by(dimension) %>%
      slice(1:ntop) %>%
      mutate(name = paste(dimension, "pos")) %>%
      ungroup()
    sig.down <- res. %>%
      group_by(dimension) %>%
      slice(max(1, n() - ntop + 1L):n()) %>%
      mutate(name = paste(dimension, "neg")) %>%
      ungroup()
    sig <- sig.up %>%
      bind_rows(sig.down) %>%
      arrange(dimension, desc(score))
  } else {
    sig <- res. %>%
      group_by(dimension) %>%
      # head(ntop) %>%
      slice(1:ntop) %>%
      mutate(name = paste(dimension, "unsigned")) %>%
      ungroup()
  }

  sig <- mutate(sig, collection = collection_name)
  sig <- select(sig, collection, name, feature_id, everything())
  # sig <- as_facile_frame(sig, fds(x))

  out <- list(
    result = sig,
    params = list(pcs = dims, ntop = ntop))

  class(out) <- sub("Ranks$", "Signature", class(x))
  out
}

#' @noRd
#' @export
signature.FacilePcaAnalysisResult <- function(x, type = "features",
                                              signed = TRUE,
                                              dims = NULL, ntop = 20, ...) {
  signature(ranks(x, type = type, signed = signed, dims = dims, ...),
            ntop = ntop, ...)
}

# Facile API ===================================================================

#' @noRd
#' @export
samples.FacilePcaAnalysisResult <- function(x, ...) {
  x[["samples"]]
}

# Utility Functions ============================================================

#' Helper function that creates a long/tidy table of statistics over the
#' features of a PCA (loadings, weights, ranks)
#' @noRd
.fpca.feature_statistics <- function(x, ...) {
  assert_class(x, "FacilePcaAnalysisResult")

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
print.FacilePcaAnalysisResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
format.FacilePcaAnalysisResult <- function(x, ...) {
  n.features <- param(x, "ntop")
  pcv <- x$percent_var * 100
  pcvs <- paste(names(pcv), sprintf("%.2f%%", pcv), sep = ": ")
  pcvu <- paste(head(pcvs, 5), collapse = "\n  ")

  out <- paste(
    "===========================================================\n",
    sprintf("FacilePcaAnalysisResult\n"),
    "-----------------------------------------------------------\n",
    "Assay used: ", param(x, "assay_name"), "\n",
    "Number of features used: ", n.features, "\n",
    "Number of PCs: ", length(pcv), "\n",
    "Variance explained:\n  ", pcvu, "\n",
    "===========================================================\n",
    sep = "")
  out
}

