#' Performs Gene Set Enrichment Analyses
#'
#' **For now** only GSEA methods process a pre-ranked feature set work, like
#' `"cameraPR"` or `"fsea"`.
#'
#' @section Updates required to multiGSEA:
#' I need to update multiGSEA to take an input data.frame of differential
#' expression statistics to work with goseq as well such that it doesn't have to
#' run the differential expression stuff again.
#'
#' Once this is implemented, the a call to fsea will also work on an Anova
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
#' ttest.gsea <- fsea(ttest.res, gdb)
#'
#' # GSEA from ANOVA result ----------------------------------------------------
#' # Not yet implemented, requires small update to multiGSEA
#'
#' # GSEA over loadings on a Principal Component -------------------------------
#' # Not yet implemented, need to get a signed weight out of eigenWeightedMean
#' # right now we just have weights. Or, fully extracing the biplot code, I
#' # think the weight should be the length on the PC of choice, and the sign is
#' # the same.
fsea <- function(x, gdb, methods, ...) {
  UseMethod("fsea", x)
}

#' @noRd
#' @export
fsea.FacileTtestAnalysisResult <- function(x, gdb, methods = "cameraPR",
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
  mg <- multiGSEA(gdb, ranked, method = methods, ...)

  # Hack the mg result until we can send in a data.frame with ranking statistics
  # and metadata
  mg@logFC <- local({
    lfc <- mg@logFC
    xref <- match(lfc[["featureId"]], ranks.[["feature_id"]])
    ranks. <- ranks.[xref,]
    stopifnot(all.equal(lfc[["featureId"]], ranks.[["feature_id"]]))
    xfer <- c("logFC", "t", "F", "B", "AveExpr", "pval", "padj",
              "meta", "symbol")
    xfer <- intersect(xfer, colnames(ranks.))
    for (cname in xfer) lfc[[cname]] <- ranks.[[cname]]
    lfc
  })

  out <- list(
    result = mg,
    params = list(methods = methods, min_logFC = min_logFC, max_padj = max_padj,
                  rank_by = rank_by))
  class(out) <- c("FacileTtestSeaAnalysisResult",
                  "FacileSeaAnalysisResult",
                  "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
fsea.FacileAnovaAnalysisResult <- function(x, gdb, methods = "goseq",
                                           max_padj = 0.10, ...) {

}

#' @noRd
#' @export
fsea.FacilePcaAnalysisResult <- function(x, methods = "cameraPR", pc = x[["pcs"]][1L],
                                         ...) {
  ranks. <- ranks(x, ...)
}

#' @section Accessing Results:
#' A call to `result(fgea.res)` will return the `multiGSEA::MultiGSEAResult`
#' object. If the user specifies the name of a GSEA method used in the analysis,
#' then the summarized results from that method will be returned by passing
#' through to the `multiGSEA::result`function, ie. `r1` and `r2` in the code
#' below are equivalent:
#'
#' ```r
#' r1 <- result(fsea.res, name = "cameraPR")
#' r2 <- result(fsea.res) %>% multiGSEA::result("cameraPR")
#' ```
#'
#' @rdname fgsea
#' @export
result.FacileGseaResult <- function(x, name = "result", ...) {
  mgres <- x[["result"]]
  if (name == "result") {
    out <- mgres
  } else {
    out <- multiGSEA::result(mgres, name, ...)
  }
  out
}
