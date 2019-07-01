#' Performs Gene Set Enrichment Analyses
#'
#' **For now** only GSEA methods process a pre-ranked feature set work, like
#' `"cameraPR"` or `"fgsea"`.
#'
#' The default method run is simply `"cameraPR"`.
#'
#' @section Updates required to multiGSEA:
#' I need to update multiGSEA to take an input data.frame of differential
#' expression statistics to work with goseq as well such that it doesn't have to
#' run the differential expression stuff again.
#'
#' Once this is implemented, the a call to `ffsea` will also work on an Anova
#' result, as well.
#'
#' @export
#' @importFrom multiGSEA multiGSEA
#'
#' @param x A `FacileAnalysisResult` object
#' @param gdb A `multiGSEA::GeneSetDb` object
#' @param methods the GSEA methods to use on `x`.
#' @return A FacileGSEAResult object, which includes a MultiGSEAResult object
#'   as it's `result()`.
#'
#' @examples
#' gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
#'
#' # GSEA from t-test result ---------------------------------------------------
#' ttest.res <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor", denom = "normal", fixed = "sex") %>%
#'   fdge(method = "voom")
#' ttest.gsea <- ffsea(ttest.res, gdb, methods = c("cameraPR", "fgsea"))
#'
#' mgsea.result <- result(ttest.gsea)
#' camera.stats <- tidy(ttest.gsea, "cameraPR")
#' fgsea.stats <- tidy(ttest.gsea, "fgsea")
#'
#' # GSEA from ANOVA result ----------------------------------------------------
#' # Not yet implemented, requires small update to multiGSEA
#'
#' # GSEA over loadings on a Principal Component -------------------------------
#' pca.crc <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' pca1.gsea <- ffsea(pca.crc, gdb, pc = 1)
#'
#' # Not yet implemented, need to get a signed weight out of eigenWeightedMean
#' # right now we just have weights. Or, fully extracing the biplot code, I
#' # think the weight should be the length on the PC of choice, and the sign is
#' # the same.
ffsea <- function(x, ...) {
  UseMethod("ffsea", x)
}

#' @noRd
#' @export
ffsea.FacileTtestAnalysisResult <- function(x, gdb, methods = "cameraPR",
                                            min_logFC = 1, max_padj = 0.10,
                                            rank_by = "logFC", signed = TRUE,
                                            ...) {
  assert_class(gdb, "GeneSetDb")
  assert_subset(methods, c("cameraPR", "fgsea", "geneSetTest"))

  ranks. <- result(ranks(x, signed = signed, ...))
  assert_choice(rank_by, colnames(ranks.))
  assert_numeric(ranks.[[rank_by]], any.missing = FALSE)

  # I can't get `arrange(ranks, desc(!!rank_by))` to work
  ranks. <- arrange_at(ranks., rank_by, desc)
  ranked <- setNames(ranks.[[rank_by]], ranks.[["feature_id"]])

  xmeta. <- select(ranks., feature_id, symbol, meta, logFC, t, B, AveExpr,
                   pval, padj)

  mg <- multiGSEA(gdb, ranked, method = methods, ..., xmeta. = xmeta.)

  out <- list(
    result = mg,
    params = list(methods = methods, min_logFC = min_logFC, max_padj = max_padj,
                  rank_by = rank_by))
  class(out) <- c("FacileTtestFseaAnalysisResult",
                  "FacileFseaAnalysisResult",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
ffsea.FacileAnovaAnalysisResult <- function(x, gdb, methods = "goseq",
                                            max_padj = 0.10, ...) {
  stop("Not yet implemented")
}

#' feature set enrichment analysis only works on one PC at a time.
#'
#' @noRd
#' @export
ffsea.FacilePcaAnalysisResult <- function(x, gdb, pc = 1,
                                          signed = TRUE, methods = "cameraPR",
                                          ...) {
  fds. <- assert_facile_data_store(fds(x))
  aname. <- assert_choice(param(x, "assay_name"), assay_names(fds.))

  messages <- character()
  warnings <- character()
  errors <- character()

  clazz <- "FacilePcaFseaAnalysisResult"
  classes <- c("FacileFseaAnalysisResult", "FacileAnalysisResult")

  out <- list(
    params = list(pc = pc, signed = signed, methods = methods, pca_result = x),
    fds = fds.)

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, classes)
    return(out)
  })


  rank.column <- if (signed) "score" else "weight"
  pc.ranks <- result(ranks(x, pcs = pc, signed = signed, ...))
  assert_choice(rank.column, colnames(pc.ranks))
  vals <- assert_numeric(pc.ranks[[rank.column]])
  names(vals) <- pc.ranks[["feature_id"]]

  xmeta. <- fds. %>%
    assay_feature_info(aname.) %>%
    transmute(feature_id, symbol = name, biotype = meta)

  mgres <- multiGSEA(gdb, vals, methods = methods, ..., xmeta. = xmeta.)
  out[["result"]] <- mgres
  out
}

#' @section Accessing Results:
#' We are in a bit of a schizophrenic state right now, where `tidy()` is
#' being the de-facto way to answer "tidy" like results (instead of result()).
#'
#' This is not to say that `result()` can't also return something that's
#' "tidy", but in this case, result(ffsea.result) will return the
#' MultiGSEAResult object itself, and `tidy(ffsea.result)` will dispatch
#' to [multiGSEA::result()] to fetch the gsea statistcs for the method
#' requested.
#'
#' ```
#' mgres <- result(ffsea.res) # return the MultiGSEAResult object
#' camera.stats <- tidy(ffsea.res, name = "cameraPR")
#' ```
#'
#' @rdname ffsea
#' @export
#' @importFrom multiGSEA resultNames
result.FacileFseaAnalysisResult <- function(x, name = "object", ...) {
  mgres <- x[["result"]]
  if (name == "object") return(mgres)

  name. <- assert_choice(name, param(x, "methods"))
  out <- as.tbl(result(mgres, name.))
  out
}

#' @noRd
#' @export
tidy.FacileFseaAnalysisResult <- function(x, name = param(x, "methods")[1L],
                                          ...) {
  mgres <- x[["result"]]
  name. <- assert_choice(name, param(x, "methods"))
  out <- as.tbl(result(mgres, name.))
  out
}

# Ranks and Signatures =========================================================

# Not sure if ranks() of a gsea analysis result should be genesets, but
# here we go.

#' @noRd
#' @export
ranks.FacileFseaAnalysisResult <- function(x, name = param(x, "methods")[1L],
                                           signed = FALSE, ...) {
  name. <- assert_choice(name, param(x, "methods"))
  rnks <- tidy(x, name.)
  if (signed) {
    rnks <- arrange(rnks, des(mean.logFC.trim))
  } else {
    rnks <- arrange(rnks, pval)
  }

  rnks <- select(rnks, collection, name, n, pval, padj,
                 padj.by.collection, everything())

  out <- list(
    result = rnks,
    params = list(name = name, signed = signed))
  # todo: need to add feature_type
  clazz <- "FacileFeatureSetRanks%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c(s, ""))
  class(out) <- classes
  out
}

#' @noRd
#' @export
result.FacileFeatureSetRanks <- function(x, name = "result") {
  x[["result"]]
}

#' @noRd
#' @export
tidy.FacileFeatureSetRanks <- function(x, name = "result") {
  x[["result"]]
}
