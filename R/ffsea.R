#' Performs Feature (Gene) Set Enrichment Analyses
#'
#' @description
#' Currently we support running a number of feature (gene) set enrichment
#' analyses downstream of a *some other* `FacileAnalysisResult` (ie.,
#' `[fdge() | fpca()] %>% ffsea()`), or over an arbitrary data.frame of
#' feature-level statistics.
#'
#' Please refer to the examples here, as well as theh *"Feature Set Analysis"*
#' section of the vignette for more information.
#'
#' @details:
#' When running `ffsea` over a `FacileAnalysisResult`, the types of methods
#' that can be run, and their configuration are preconfigured with reaonable
#' defaults.
#'
#' When providing a generic `data.frame` of feature-level statistics to run
#' enrichment tests over, the user has to specificy a few more parmaters.
#' If running a "pre-ranked" test, (`method %in% c("cameraPR", "fgsea")`) the
#' name of a numeric column in `x` must be specified to rank the features by,
#' and the order by which to do that using the `rank_by` and `rank_order`
#' arguments, respectively.
#'
#' If running an overrepresentation analysis-style test (`method = "ora"`),
#' the user must specifcy the name of a logical column that indicates
#' (when `TRUE`) that a feature should be included for enrichment testing. The
#' user can optionally specify a `group_by` column, like `"direction"`, that
#' will be used to split the selected features into groups to perform more
#' specific enrichmen tests. This allows you enrichment tests to be run
#' separately for `"up"` and `"down"` regulated genes separately, for example.
#'
#' Lastly, the user can provide the name of another numeric column in `x` with
#' `biased_by` which can be used to account for bias in the enrichment tests,
#' such as gene length, GC content, etc.
#'
#' Gene sets must be supplied as a [multiGSEA::GeneSetDb()] object.
#'
#' @section GSEA Methods:
#' Currently, only the following GSEA methods are supported:
#'
#' * `"cameraPR"`: Delegates to [limma::cameraPR()] to perform a competitive
#'   gene set test based on feature ranks imposed downstream of an analysis
#' * `"fgsea"`: Delegates to [fgsea::fgsea()] to perform another version of
#'   a competitive gene set test based on ranks.
#' * `"ora"`: Performs an overrepresentation analysis test. The user
#'   must specify the name of `logical` column (`select_by`) from the input
#'   which is used to indicate the features that are selected for enrichment
#'   analysis. The user can optionally provide the name of a `numeric` column
#'   (`biased_by`) and `character` column (`group_by`), which will adjust the
#'   enrichment test for a covariate that may induce a bias in the DGE results,
#'   and also run follow up enrichment tests based by differnt groups of
#'   features (`group_by`). For example, the result table might have a
#'   `"direction"` column, which specifies the direciton of differential
#'   expression (`"up"`, or `"down"`). In this case, enrichment tests will be
#'   run over *all* features together, and then independantly for the ones that
#'   are `"up"`, and `"down"`.
#'
#' @section GSEA Statistics:
#' The geneset level statistics can be extracted from the
#' `FacileFseaAnalysisResult` on a per-method basis usig the `tidy()` function.
#' For instance, if `ffsea()` was called with
#' `fres <- ffsea(..., methods = c("cameraPR", "ora")`, the `"cameraPR"`
#' results can be extracted via `tidy(fres, "cameraPR")`
#'
#' @section Development Notes:
#' This functionality delegates to multiGSEA to do all of the work. The
#' multiGSEA interface is undergoing a bit of refactoring in order to better
#' support a table of feature statistcs as input (for preranked and enrichment
#' tests), so the `"methods"` supported via `ffsea()` are limited to a subset
#' of the ones wrapped by multiGSEA, as enumerated below.

#' @export
#' @importFrom multiGSEA multiGSEA
#' @seealso https://github.com/lianos/multiGSEA
#' @aliases gsea GSEA fsea FSEA
#'
#' @param x A `FacileAnalysisResult` object, or a data.frame with feature-level
#'   statistics, minimally with a `"feature_id"` column as well as one or more
#'   `numeric` columns to rank features on.
#' @param gdb A [multiGSEA::GeneSetDb()] object
#' @param methods the GSEA methods to use on `x`.
#' @return A FacileFseaAnalysisResult object, which includes a MultiGSEAResult
#'   object as it's `result()`. The geneset level statistics for each of the
#'   methods that were run are available via `tidy(ffsea.res, "<method_name>")`.
#'
#' @examples
#' gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
#' efds <- FacileData::exampleFacileDataSet()
#'
#' # GSEA from t-test result ---------------------------------------------------
#' ttest.res <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
#'           batch = "sex") %>%
#'   fdge(method = "voom")
#'
#' ttest.gsea <- ffsea(ttest.res, gdb, methods = c("cameraPR", "ora"),
#'                     biased_by = "effective_length")
#' if (interactive()) {
#'   ttest.gsea <- ffseaGadget(ttest.res, gdb)
#'   shine(ttest.gsea)
#' }
#'
#' camera.stats <- tidy(ttest.gsea, "cameraPR")
#' ora.stats <- tidy(ttest.gsea, "ora")
#'
#' # GSEA from ANOVA result ----------------------------------------------------
#' stage.anova <- efds %>%
#'   FacileData::filter_samples(indication == "BLCA") %>%
#'   flm_def(covariate = "stage", batch = "sex") %>%
#'   fdge(method = "voom")
#' anova.gsea <- ffsea(stage.anova, gdb)
#' if (interactive()) {
#'  shine(anova.gsea)
#'  # We can generate the same GSEA result like so
#'  anova.gsea2 <- ffseaGadget(stage.anova, gdb = gdb)
#' }
#' # GSEA over loadings on a Principal Component -------------------------------
#' pca.crc <- efds %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca()
#' pca1.gsea <- ffsea(pca.crc, gdb, dim = 1)
ffsea <- function(x, gdb, methods = "cameraPR", ...) {
  UseMethod("ffsea", x)
}

#' TODO: ffsea.FacileTtestDGEModelDefinition has all the info required to run
#' GSEA methods that require more than just pre-ranked (or enriched) listing
#' of genes.
#'
#' @noRd
ffsea.FacileTtestDGEModelDefinition <- function(x, gdb, methods = "cameraPR",
                                                ...) {
  stop("Run fdge(x) %>% ffsea() for now")
}

#' @section Generic Set Enrichment Analysis from a table of statistics:
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
#'   over represenatation tests. Rows where x[[select_by]] is `TRUE` are
#'   included for enrichment analysis
#' @param rank_order the direction to arrange values in `rank_by`. By default
#'   (`rank_by = "asc"`), which arranges `x[[rank_by]]` in ascending order.
#'   Specifying `rank_by = "desc"` will rank `x` by `rank_by` in descending
#'   order. If `"rankded"`, then we assume that the data.frame is already
#'   ranked as desired.
ffsea.data.frame <- function(x, gdb, methods = "cameraPR",
                             rank_by = NULL, select_by = NULL,
                             rank_order = "descending", group_by = NULL,
                             biased_by = NULL, ...,
                             feature.bias = "thebuckstopshere",
                             xmeta. = "thebuckstopshere",
                             groups = "thebuckstopshere") {
  if (is.factor(x[["feature_id"]])) {
    x[["feature_id"]] <- as.character(x[["feature_id"]])
  }
  assert_character(x[["feature_id"]])
  assert_class(gdb, "GeneSetDb")

  methods. <- local({
    assert_character(methods, min.len = 1L)
    m <- filter(.ffsea_methods(), method %in% methods)
    assert_set_equal(m[["method"]], methods)
    m
  })

  types <- unique(methods.[["type"]])

  xx <- x

  if ("ranks" %in% types) {
    assert_choice(rank_by, colnames(x))
    assert_numeric(x[[rank_by]])
    if (missing(rank_order)) {
      if (rank_by %in% c("pval", "padj")) {
        rank_order <- "ascending"
      } else {
        rank_order <- "descending"
      }
    }
    rank_order <- match.arg(rank_order, c("ascending", "descending", "ranked"))
    arrange.fn <- if (rank_order == "descending") dplyr::desc else list()
    if (rank_order != "ranked") {
      xx <- arrange_at(xx, rank_by, arrange.fn)
    }
  }

  if ("enrichment" %in% types) {
    assert_choice(select_by, colnames(x))
    assert_logical(xx[[select_by]])
    if (select_by != "significant") {
      # currently multiGSEA() enrichment methods like goseq or hyperG require
      # a "significant" and "direction" column
      xx[["significant"]] <- xx[[select_by]]
    }
    if (!is.null(group_by)) {
      assert_choice(group_by, colnames(x))
      assert_character(xx[[group_by]])
    }
    if (!is.null(biased_by)) {
      assert_choice(biased_by, colnames(x))
      assert_numeric(xx[[biased_by]])
    }
  }

  if ("ranks" %in% types) {
    input <- xx[[rank_by]]
  } else {
    input <- as.integer(xx[[select_by]])
  }
  names(input) <- xx[["feature_id"]]

  mg <- multiGSEA(gdb, input, methods = methods, feature.bias = biased_by,
                  groups = group_by, xmeta. = xx, ...)

  out <- list(
    result = mg,
    params = list(methods = methods,
                  rank_by = rank_by, select_by = select_by,
                  rank_order = rank_order,
                  biased_by = biased_by,
                  x = x, gdb = gdb))
    fds = suppressWarnings(fds(x))
  class(out) <- c("FacileFseaAnalysisResult", "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
#' @importFrom multiGSEA multiGSEA
ffsea.FacileTtestAnalysisResult <- function(x, gdb,
                                            methods = c("cameraPR", "ora"),
                                            min_logFC = param(x, "treat_lfc"),
                                            max_padj = 0.10,
                                            rank_by = "logFC",
                                            signed = TRUE,
                                            biased_by = NULL, ...,
                                            rank_order = "ranked",
                                            group_by = "direction",
                                            select_by = "significant") {
  if (assert_string(rank_order) != "ranked") {
    warning("non-default value used for `rank_order`")
  }
  if (assert_string(group_by) != "direction") {
    warning("non-default value used for `group_by`")
  }
  if (assert_string(select_by) != "significant") {
    warning("non-default value used for `select_by`")
  }
  assert_class(gdb, "GeneSetDb")
  all.methods <- ffsea_methods(x)
  assert_subset(methods, all.methods[["method"]], empty.ok = FALSE)
  fds. <- assert_facile_data_store(fds(x))
  if (is.null(min_logFC)) min_logFC <- 0
  assert_number(min_logFC, lower = 0, finite = TRUE)

  ranks. <- tidy(ranks(x, signed = signed, rank_by = rank_by, ...))

  ranks. <- mutate(ranks.,
                   significant = padj <= max_padj, abs(logFC) >= min_logFC,
                   direction = ifelse(logFC > 0, "up", "down"))

  take.cols <- c(
    "significant", "direction",
    "symbol", "meta", "logFC", "t", "B",
    "AveExpr", "pval", "padj", "CI.L", "CI.R", "effective_length")
  take.cols <- intersect(take.cols, colnames(ranks.))

  input <- select(ranks., feature_id, {{take.cols}})

  messages <- character()
  warnings <- character()
  errors <- character()

  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(
      c("FacileTtestFseaAnalysisResult", "FacileDgeFseaAnalysisResult"),
      class(out))
    return(out)
  })

  out <- ffsea(input, gdb, methods = methods,
               rank_by = rank_by, rank_order = "ranked",
               select_by = "significant", group_by = "direction",
               biased_by = biased_by, ...)

  out[["params"]][["xdf"]] <- out[["params"]][["x"]]
  out[["params"]][["x"]] <- x
  out[["params"]][["min_logFC"]] <- min_logFC
  out[["params"]][["max_padj"]] <- max_padj
  out[["params"]][["signed"]] <- signed
  out[["fds"]] <- fds(x)
  out
}

#' @rdname fdge
#' @export
#' @section Gene Set Enrichment Analysis:
#' There are a few ways you may consider running a gene set analysis over an
#' interaction analysis.
#'
#' 1. On the statistics of the interaction itself; or
#' 2. On the statistics of the different "quadrants" the features are binned
#'    into that are found in the `sigclass` columns of the
#'    tidied result table, ie. `tidy.FacileTtestComparisonAnalysisResult()`; or
#' 3. Both.
#'
#' Note that an analysis on (2) only lends itself to an overrepresentation
#' analysis, ie. `methods = "ora"`.
ffsea.FacileTtestComparisonAnalysisResult <- function(
    x, gdb, methods = c("cameraPR", "ora"),
    type = c("interaction", "quadrants"),
    min_logFC = param(x, "treat_lfc"),
    min_logFC_x = min_logFC, min_logFC_y = min_logFC,
    max_padj = 0.10, max_padj_x = max_padj, max_padj_y = max_padj,
    rank_by = "logFC", signed = TRUE, biased_by = NULL, ...,
    rank_order = "ranked", group_by = "direction", select_by = "significant") {
  assert_class(gdb, "GeneSetDb")
  type <- match.arg(type)
  if (type == "interaction") {
    # NextMethod should call ffsea.FacileTtestAnalysisResult
    # out <- NextMethod(override = parameters, as_you = like)
    out <- NextMethod()
    class(out) <- c("FacileTtestComparisonInteractionFseaAnalysisResult")
  } else {
    out <- .ffsea.ttestcomp.quadrants()
    class(out) <- c("FacileTtestComparisonQuadrantFseaAnalysisResult")
  }

  out
}

#' @noRd
.ffsea.ttest.comp.interaction <- function(){
  params <- list(
    x = x, gdb = gdb, methods = methods, type = type,
    min_logFC = min_logFC, min_logFC_x = min_logFC_x, min_logFC_y = min_logFC_y,
    max_padj = max_padj, max_padj_x = max_padj_x, max_padj_y = max_padj_y,
    rank_by = rank_by, signed = signed, biased_by = biased_by,
    rank_order = rank_order, group_by = group_by, select_by = select_by)

  out <- list(
    params = params)
}

#' @noRd
.ffsea.ttest.comp.quadrants <- function() {
  params <- list(
    x = x, gdb = gdb, methods = methods, type = type,
    min_logFC = min_logFC, min_logFC_x = min_logFC_x, min_logFC_y = min_logFC_y,
    max_padj = max_padj, max_padj_x = max_padj_x, max_padj_y = max_padj_y,
    rank_by = rank_by, signed = signed, biased_by = biased_by,
    rank_order = rank_order, group_by = group_by, select_by = select_by)

  out <- list(
    params = params)
}

#' @noRd
#' @export
ffsea.FacileAnovaAnalysisResult <- function(x, gdb, methods = "ora",
                                            max_padj = 0.10, biased_by = NULL,
                                            ...,
                                            rank_by = "F",
                                            select_by = "significant",
                                            rank_order = "ranked") {
  if (assert_string(rank_by) != "F") {
    warning("non-default value used for `rank_order`")
  }
  if (assert_string(select_by) != "significant") {
    warning("non-default value used for `select_by`")
  }
  if (assert_string(rank_order) != "ranked") {
    warning("non-default value used for `group_by`")
  }

  all.methods <- ffsea_methods(x)
  assert_subset(methods, all.methods[["method"]], empty.ok = FALSE)
  fds. <- assert_facile_data_store(fds(x))

  ranks. <- tidy(x)
  ranks. <- mutate(ranks., significant = padj <= max_padj)

  out <- ffsea(ranks., gdb, methods = methods,
               select_by = select_by, rank_by = rank_by,
               rank_order = rank_order, biased_by = biased_by, ...)
  out[["params"]][["xdf"]] <- out[["params"]][["x"]]
  out[["params"]][["x"]] <- x
  out[["params"]][["max_padj"]] <- max_padj
  out[["fds"]] <- fds(x)

  class(out) <- c(
    c("FacileAnovaFseaAnalysisResult", "FacileDgeFseaAnalysisResult"),
    class(out))
  out
}

#' feature set enrichment analysis only works on one PC at a time.
#'
#' @noRd
#' @export
#' @importFrom multiGSEA multiGSEA
ffsea.FacilePcaAnalysisResult <- function(x, gdb, methods = "cameraPR", dim = 1,
                                          signed = TRUE, ...) {
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
    # class(out) <- c(clazz, classes)
    class(out) <- c(clazz, class(out))
    return(out)
  })

  rank.column <- if (signed) "score" else "weight"
  pc.ranks <- tidy(ranks(x, dims = dim, signed = signed, ...))

  out <- ffsea(pc.ranks, gdb, methods = methods, rank_by = rank.column,
               rank_order = "descending", ...)

  out[["params"]][["xdf"]] <- out[["params"]][["x"]]
  out[["params"]][["x"]] <- x
  out[["params"]][["dim"]] <- dim
  out[["params"]][["signed"]] <- signed
  out[["fds"]] <- fds(x)
  out
}

# Methods and Accessors ========================================================

#' @noRd
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
  out <- as_tibble(result(mgres, name.))
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
  # name. <- assert_choice(name, param(x, "methods"))
  name. <- assert_choice(name, multiGSEA::resultNames(mgres))
  out <- as_tibble(result(mgres, name.))
  select(out, collection, name, pval, padj, padj.by.collection, everything())
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
