#' Performs Feature (Gene) Set Enrichment Analyses
#'
#' **For now** only GSEA methods process a pre-ranked feature set work, like
#' `"cameraPR"` and `"fgsea"`.
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
#' efds <- FacileData::exampleFacileDataSet()
#'
#' # GSEA from t-test result ---------------------------------------------------
#' ttest.res <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor", denom = "normal", fixed = "sex") %>%
#'   fdge(method = "voom")
#' ttest.gsea <- ffsea(ttest.res, gdb, methods = c("cameraPR", "fgsea"))
#' if (interactive()) {
#'   shine(ttest.gsea)
#' }
#'
#' mgsea.result <- result(ttest.gsea)
#' camera.stats <- tidy(ttest.gsea, "cameraPR")
#' fgsea.stats <- tidy(ttest.gsea, "fgsea")
#'
#' # GSEA from ANOVA result ----------------------------------------------------
#' \dontrun{
#' # This requires an update in mutiGSEA (to support df inputs) to run
#' stage.anova <- efds %>%
#'   FacileData::filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "stage", fixed = "sex") %>%
#'   fdge(method = "voom")
#' anova.gsea <- ffsea(stage.anova)
#' if (interactive()) {
#'  shine(anova.gsea)
#' }
#' }
#'
#' # GSEA over loadings on a Principal Component -------------------------------
#' pca.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' pca1.gsea <- ffsea(pca.crc, gdb, dim = 1)
#'
#' # Not yet implemented, need to get a signed weight out of eigenWeightedMean
#' # right now we just have weights. Or, fully extracing the biplot code, I
#' # think the weight should be the length on the PC of choice, and the sign is
#' # the same.
ffsea <- function(x, ...) {
  UseMethod("ffsea", x)
}

#' TODO: ffsea.FacileTtestDGEModelDefinition has all the info required to run
#' GSEA methods that require more than just pre-ranked (or enriched) listing
#' of genes.
#'
#' @noRd
ffsea.FacileTtestDGEModelDefinition <- function(x, gdb,
                                                methods = c("camera", "roast"),
                                                ...) {
}

#' @section Generic Set Enrichment Analysis from a table of statistics:
#' TODO: Write generic ffsea over a data.frame of feature statistics. This can
#'   drive pre-ranked based GSEA, by identifying the column in `x` to `rank_by`,
#'   as well as "enrichment" based methods, by identifying an indicator column
#'   `select_by` in `x` flags features as "interesting" or not.
#'
#' TODO: Implement the universal ffsea.data.frame
#'
#' @noRd
#' @export
#' @method ffsea data.frame
#'
#' @param x a data.frame, rows are features, columns are metadata or statistics
#'   over the features
#' @param rank_by the name of a numeric column in `x` to use to arrange the
#'   features by ranks
#' @param select_by a logical column in `x` used to select features for
#'   enrichmen tests. Rows where x[[select_by]] is `TRUE` are included for
#'   enrichment analysis
#' @param rank_order the direction to arrange values in `rank_by`. By default
#'   (`rank_by = "asc"`), which arranges `x[[rank_by]]` in ascending order.
#'   Specifying `rank_by = "desc"` will rank `x` by `rank_by` in descending
#'   order.
ffsea.data.frame <- function(x, rank_by, select_by, methods,
                             rank_order = c("asc", "desc"), ...) {

}

#' @noRd
#' @export
#' @importFrom multiGSEA multiGSEA
ffsea.FacileTtestAnalysisResult <- function(x, gdb, methods = "cameraPR",
                                            min_logFC = 1, max_padj = 0.10,
                                            rank_by = "logFC", signed = TRUE,
                                            ...) {
  assert_class(gdb, "GeneSetDb")
  assert_subset(methods, c("cameraPR", "fgsea", "geneSetTest"))
  fds. <- assert_facile_data_store(fds(x))

  ranks. <- result(ranks(x, signed = signed, ...))
  assert_choice(rank_by, colnames(ranks.))
  assert_numeric(ranks.[[rank_by]], any.missing = FALSE)

  # I can't get `arrange(ranks, desc(!!rank_by))` to work
  ranks. <- arrange_at(ranks., rank_by, desc)
  ranked <- setNames(ranks.[[rank_by]], ranks.[["feature_id"]])

  messages <- character()
  warnings <- character()
  errors <- character()

  clazz <- c("FacileTtestFseaAnalysisResult", "FacileDgeFseaAnalysisResult")
  classes <- c("FacileFseaAnalysisResult", "FacileAnalysisResult")

  out <- list(
    result = NULL,
    params = list(methods = methods, min_logFC = min_logFC,
                  max_padj = max_padj, rank_by = rank_by, x = x),
    fds = fds.)

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, classes)
    return(out)
  })

  take.cols <- c("symbol", "meta", "logFC", "t", "B", "AveExpr", "pval", "padj",
                 "CI.L", "CI.R")
  take.cols <- intersect(take.cols, colnames(ranks.))
  xmeta. <- select(ranks., feature_id, {{take.cols}})

  out[["result"]] <- multiGSEA(gdb, ranked, method = methods, ...,
                               xmeta. = xmeta.)
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
#' @importFrom multiGSEA multiGSEA
ffsea.FacilePcaAnalysisResult <- function(x, gdb, dim = 1,
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
    params = list(dim = dim, signed = signed, methods = methods, x = x),
    fds = fds.)

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, classes)
    return(out)
  })

  rank.column <- if (signed) "score" else "weight"
  pc.ranks <- result(ranks(x, dims = dim, signed = signed, ...))
  assert_choice(rank.column, colnames(pc.ranks))
  vals <- assert_numeric(pc.ranks[[rank.column]])
  names(vals) <- pc.ranks[["feature_id"]]

  # xmeta. <- fds. %>%
  #   features(aname.) %>%
  #   transmute(feature_id, symbol = name, score, weight, biotype = meta)

  xmeta. <- pc.ranks
  out[["result"]] <- multiGSEA(gdb, vals, methods = methods, ...,
                               xmeta. = xmeta.)
  out
}

# Methods and Accessors ========================================================

#' @section Feature Set Enrichment Analysis:
#' What are the features of a feature-set enrichment analysis ([ffsea()])?
#' Aren't they the gene sets, and not the individual genes themselves?
#' There is crappy support for this, for now.
#'
#' It is a meta-something type of thing. The genesets are the features, but
#' they are also made up of their own features. We most often think of genesets
#' as consisting of genes, but perhaps we can imagine a feature set that
#' consists of motifs ... or sometihng.
#'
#' @rdname features
#' @export
#' @importFrom multiGSEA encode_gskey
features.FacileFseaAnalysisResult <- function(x, ...) {
  warning("The feature_id,feature_type feature representation for fsea is a ",
          "bit loose, refer to the 'Feature Set Enrichment Analyais' section ",
          "of ?features")?
  stat.table <- tidy(x)
  stat.table[["feature_id"]] <- encode_gskey(stat.table)
  stat.table[["feature_type"]] <- "feature_set"
  select(stat.table, collection, name, feature_id, feature_type, everything())
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
  if (name == "object") {
    return(mgres)
  }

  name. <- assert_choice(name, param(x, "methods"))
  out <- as.tbl(result(mgres, name.))
  out
}

#' @noRd
#' @export
initialized.FacileFseaAnalysisResult <- function(x, ...) {
  is(result(x), "MultiGSEAResult")
}

#' @noRd
#' @export
samples.FacileFseaAnalysisResult <- function(x, ...) {
  x.parent <- param(x, "x")
  samples(x.parent)
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

# Printing =====================================================================

#' @noRd
#' @export
print.FacileFseaAnalysisResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
#' @importFrom multiGSEA resultNames tabulateResults
format.FacileFseaAnalysisResult <- function(x, max_padj = 0.20, ...) {
  mgres <- result(x)
  gsea.res.table <- tabulateResults(mgres, max.p = max_padj)
  source.type <- class(param(x, "x"))[1L]

  if (source.type == "FacilePcaAnalysisResult") {
    source.type <- sprintf("%s [PC: %s]", source.type,
                           as.character(param(x, "dim")))
  }

  msg <- paste(
    paste(rep("=", 80), collapse = ""), "\n",
    sprintf("FacileFseaAnalysisResult (from a %s)\n", source.type),
    paste(rep("-", 80), collapse = ""), "\n",
    paste(tibble:::format.tbl(gsea.res.table)[-1], collapse = "\n"), "\n",
    paste(rep("=", 80), collapse = ""), "\n",
    sep = "")
  msg
}
