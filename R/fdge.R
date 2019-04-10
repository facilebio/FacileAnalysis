#' Peform a differential expression analysis.
#'
#' Use [fdge_model_def()] to define the design matrix and contrast to test and
#' pass the `FacileDGEModelDefinition` object returned from that to `fdge()`
#' to run the desired differential testing framework (dictated by the `method`
#' parameter) over the data. `fdge_model_def` accepts a
#'
#' @section Differential Expression Testing Methods:
#' The appropriate statistical framework to use for differential expression
#' testing is defined by the type of data that is recorded in the assay
#' `assay_name`, ie. `assay_info(x, assay_name)$assay_type`.
#'
#' The `fdge_methods()` function returns a tibble of appropriate
#' `assay_type -> dge_method` associations. The first entry for each
#' `dge_method` is the default `method` used if one isn't provided by the
#' caller.
#'
#' The available methods are:
#'
#' * `"voom"`: For count data, uses [limma::voomWithQualityWeights()] when
#'   `with_sample_weights = TRUE`.
#' * `"edgeR-qlf"`: The edgeR quasi-likelihood method, for count data.
#' * `"limma-trend"`: Usable for log-transformed data that "looks like" it came
#'    from count data, or where there is a "trend" of the variance with the
#'    mean, uses  [lima::arrayWeights()] when `with_sample_weights = TRUE`.
#' * `"limma"`: Straightup limma, this expects log2-normal like data, with
#'    (largely) no trend of the variance to the mean worth modeling. Uses
#'    [lima::arrayWeights()] when `with_sample_weights = TRUE`
#'
#' @export
#' @importFrom multiGSEA GeneSetDb calculateIndividualLogFC logFC multiGSEA
#'
#' @param x a data source
#' @param assay_name the name of the assay that holds the measurements for test.
#'   Defaults to `default_assay(x)`.
#' @param method The differential testing framework to use over the data. The
#'   valid choices are defined by the type of assay `assay_name`is. Refer to the
#'   *Differential Expression Testing Methods* section for more details.
#' @param filter,with_sample_weights Passed into [biocbox()] to determine which
#'   features (genes) are removed from the dataset for testing, as well as if
#'   to use [limma::arrayWeights()] or [limma::voomWithQualityWeights()]
#'   (where appropriate) when testing (default is not to).
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#' samples <- FacileData::filter_samples(efds, indication == "BLCA")
#' mdef <- fdge_model_def(samples, covariate = "sample_type",
#'                        numer = "tumor", denom = "normal")
#' dge <- fdge(mdef, method = "voom")
#' dge.sig <- signature(dge)
#' if (interactive()) {
#'   viz(dge)
#'   shine(dge)
#' }
#'
#' stage.anova <- samples %>%
#'   fdge_model_def(covariate = "stage", fixed = "sex") %>%
#'   fdge(method = "voom")
#' anova.sig <- signature(stage.anova)
fdge <- function(x, ...) {
  UseMethod("fdge", x)
}

#' @export
#' @rdname fdge
fdge.FacileAnovaModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                            filter = "default",
                                            with_sample_weights = FALSE, ...) {
  res <- NextMethod(coef = x$coef)
  res
}

#' @export
fdge.FacileTtestDGEModelDefinition <- function(x, assay_name = NULL,
                                               method = NULL,
                                               filter = "default",
                                               with_sample_weights = FALSE,
                                               treat_lfc = NULL,
                                               ...) {
  res <- NextMethod(contrast = x$contrast)
  res
}

#' @export
#' @rdname fdge
fdge.FacileDGEModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                          filter = "default",
                                          with_sample_weights = FALSE,
                                          treat_lfc = NULL,
                                          ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  test_type <- assert_choice(x[["test_type"]], c("ttest", "anova"))
  .fds <- assert_facile_data_store(fds(x))

  if (test_type == "ttest") {
    testme <- x[["contrast"]]
    clazz <- "FacileTtestDGEResult"
  } else {
    testme <- x[["coef"]]
    clazz <- "FacileAnovaDGEResult"
  }

  if (is.null(assay_name)) assay_name <- default_assay(.fds)
  ainfo <- try(assay_info(.fds, assay_name), silent = TRUE)
  if (is(ainfo, "try-error")) {
    msg <- glue("assay_name `{assay_name}` not present in FacileDataStore")
    errors <- c(errors, msg)
    assay_type <- "ERROR"
  } else {
    assay_type <- ainfo$assay_type
  }

  dge_methods <- fdge_methods(assay_type)
  if (nrow(dge_methods) == 0L) {
    msg <- glue("Differential expression method `{method}` for assay_type ",
                "`{assay_type}` not found")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0L) {
    bb <- biocbox(x, assay_name, method, dge_methods, filter,
                  with_sample_weights = with_sample_weights, ...)
    messages <- c(messages, bb[["messages"]])
    warnings <- c(warnings, bb[["warnings"]])
    errors <- c(errors, bb[["errors"]])

    method <- bb[["dge_method"]]

    if (test_number(treat_lfc)) {
      treat_lfc <- abs(treat_lfc)
      use.treat <- TRUE
    } else {
      if (!is.null(treat_lfc)) {
        warnings <- c(
          warnings,
          "Illegal parameter passed to `treat_lfc`. It is being ignored")
      }
      treat_lfc <- 1
      use.treat <- FALSE
    }

    y <- result(bb)
    des <- design(bb)
    result <- calculateIndividualLogFC(y, des, contrast = testme,
                                       use.treat = use.treat,
                                       feature.min.logFC = treat_lfc)
    # multiGSEA::calculateIndividualLogFC returns the stats table ordered by
    # featureId, let's put the features back in the order they are in y
    rownames(result) <- result[["featureId"]]
    result <- result[rownames(y),]

    axe.cols <- c("featureId", "x.idx")
    for (col in axe.cols) {
      if (col %in% names(result)) result[[col]] <- NULL
    }
    result <- as.tbl(result)
  } else {
    result <- NULL
  }

  out <- list(
    test_type = test_type,
    result = result,
    assay_name = assay_name,
    method = method,
    model_def = x,
    treat_lfc = if (use.treat) treat_lfc else NULL,
    # Standard FacileAnalysisResult things
    fds = .fds,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(out) <- c(clazz, "FacileDGEResult", "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
result.FacileDGEResult <- function(x, name, ...) {
  out <- x[["result"]]
  if (!"feature_type" %in% colnames(out)) {
    feature_type <- guess_feature_type(out[["feature_id"]], summarize = TRUE)
    out[["feature_type"]] <- feature_type[["feature_type"]]
    out <- select(feature_type, everything())
  }
  out
}

# Ranks and Signatures =========================================================

#' @export
#' @noRd
ranks.FacileTtestDGEResult <- function(x, ...) {
  ranks. <- result(x, ...)
  ranks. <- arrange(ranks., desc(logFC))
  out <- list(
    result = ranks.,
    ranking_columns = "logFC",
    ranking_order = "decreasing")
  class(out) <- c("FacileTtestDGEFeatureRankings",
                  "FacileFeatureRankings",
                  "FacileAnalysisResult")
  out
}

#' @export
#' @noRd
signature.FacileTtestDGEResult <- function(x, min_logFC = x[["treat_lfc"]],
                                           max_padj = 0.10, ntop = 20,
                                           name = NULL, collection_name = NULL,
                                           ...) {
  signature(ranks(x, ...), min_logFC = min_logFC, max_padj = max_padj,
            ntop = ntop, name = name, collection_name = collection_name, ...)
}

#' @export
#' @noRd
signature.FacileTtestDGEFeatureRankings <- function(x, min_logFC = x[["treat_lfc"]],
                                                    max_padj = 0.10,
                                                    ntop = 20,
                                                    name = NULL,
                                                    collection_name = NULL,
                                                    ...) {
  if (is.null(name)) name <- "Ttest signature"
  name. <- assert_string(name)

  if (is.null(collection_name)) collection_name <- class(x)[1L]
  assert_string(collection_name)

  if (is.null(min_logFC)) min_logFC <- 1
  up <- result(x) %>%
    filter(padj <= max_padj, logFC >= min_logFC) %>%
    mutate(collection = collection_name, name = paste(name., "up"),
           direction = "up")
  up <- head(up, ntop)

  down <- result(x) %>%
    filter(padj <= max_padj, logFC <= min_logFC) %>%
    mutate(collection = collection_name, name = paste(name., "down"),
           direction = "down")
  down <- tail(down, ntop)

  sig <- bind_rows(up, down) %>%
    select(collection, name, everything())

  out <- list(
    result = sig,
    params = list(max_padj = max_padj, min_logFC = min_logFC))
  class(out) <- sub("Rankings$", "Signature", class(x))
  out
}

#' @noRd
#' @export
ranks.FacileAnovaDGEResult <- function(x, ...) {
  ranks. <- result(x, ...)
  ranks. <- arrange(ranks., pval)
  out <- list(
    result = ranks.,
    ranking_colums = "pval",
    ranking_order = "ascending")
  class(out) <- c("FacileAnovaDGEFeatureRankings",
                  "FacileFeatureRankings",
                  "FacileAnalysisResult")
  out
}

#' @export
#' @noRd
signature.FacileAnovaDGEResult <- function(x, max_padj = 0.10, ntop = 20,
                                           name = NULL, collection_name = NULL,
                                           ...) {
  signature(ranks(x, ...), max_padj = max_padj, ntop = ntop, name = name,
            collectoin_name = collection_name, ...)
}

#' @export
#' @noRd
signature.FacileAnovaDGEFeatureRankings <- function(x, max_padj = 0.10,
                                                    ntop = 20, name = NULL,
                                                    collection_name = NULL,
                                                    ...) {
  if (is.null(name)) name <- "Anova signature"
  name. <- assert_string(name)

  if (is.null(collection_name)) collection_name <- class(x)[1L]
  assert_string(collection_name)

  sig <- result(x) %>%
    filter(padj <= max_padj) %>%
    mutate(collection = collection_name, name = name.,
           direction = "udefined") %>%
    head(ntop)

  out <- list(
    result = sig,
    params = list(ntop = ntop, max_padj = max_padj))
  class(out) <- sub("Rankings$", "Signature", class(x))
  sig
}

# Facile API ===================================================================

#' @noRd
#' @export
samples.FacileDGEResult <- function(x, ...) {
  samples(x[["model_def"]])
}

#' @section FacileDGEResult:
#' Given a FacileDGEResult, we can re-materialize the Bioconductor assay
#' container used within the differential testing pipeline used from [fdge()].
#' Currently we have limited our analysis framework to either work over DGEList
#' (edgeR) or EList (limma) containers.
#'
#' @export
#' @rdname biocbox
biocbox.FacileDGEResult <- function(x, ...) {
  res <- biocbox(x[["model_def"]], x[["assay_name"]], x[["method"]],
                 filter = result(x)[["feature_id"]], ...)
  res
}

#' @noRd
#' @export
print.FacileDGEResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacileDGEResult <- function(x, ...) {
  test_type <- if (is(x, "FacileTtestDGEResult")) "ttest" else "ANOVA"
  mdef <- x[["model_def"]]
  formula <- mdef[["design_formula"]]

  if (test_type == "ttest") {
    test <- sprintf("%s: (%s) - (%s)",
                    mdef[["covariate"]],
                    paste(mdef[["numer"]], collapse = "+"),
                    paste(mdef[["denom"]], collapse = "+"))
  } else {
    des <- mdef[["design"]]
    test <- sprintf("%s (%s)", mdef[["covariate"]],
                    paste(colnames(des)[mdef[["coef"]]], collapse = "|"))
  }
  res <- result(x)
  ntested <- nrow(res)
  nsig <- sum(!is.na(res[["padj"]]) & res[["padj"]] < 0.10)
  nsamples <- nrow(x$model_def$covariates)
  out <- paste(
    "===========================================================\n",
    sprintf("FacileDGEResult (%s)\n", test_type),
    "-----------------------------------------------------------\n",
    glue("Significant Results (FDR < 0.1): ({nsig} / {ntested})"), "\n",
    "Formula: ", formula, "\n",
    "Tested: ", test, "\n",
    "Number of samples: ", nsamples, "\n",
    "===========================================================\n",
    sep = "")
  out
}


# Helpers ======================================================================

#' A table of assay_type,dge_method combination parameters
#'
#' This table is used to match the assay_type,dge_method combination with
#' the appropriate bioc container class `bioc_class`, and default edgeR/limma
#' model fitting params.
#'
#' @export
#' @param assay_type An optional string specifying the valid dge methods for
#'   a given assay type.
#' @return A tibble of assay_type -> method and parameter associtiations. If
#'   `assay_type`is not `NULL`, this will be filtered to the associations
#'   valid only for that `assay_type`. If none are found, this will be a
#'   0-row tibble.
fdge_methods <- function(assay_type = NULL) {
  # assay_type values : rnaseq, umi, affymrna, affymirna, log2

  # This is a table of assay_type : dge_method possibilites. The first row
  # for each assay_type is the default analysis method
  assay_methods <- tribble(
    ~assay_type,   ~dge_method,         ~bioc_class,
    "rnaseq",      "voom",              "EList",
    "rnaseq",      "edgeR-qlf",         "DGEList",
    "rnaseq",      "limma-trend",       "EList",
    "umi",         "voom",              "EList",
    "umi",         "edgeR-qlf",         "DGEList",
    "umi",         "limma-trend",       "EList",
    "tpm",         "limma-trend",       "EList",
    "affymrna",    "limma",             "EList",
    "affymirna",   "limma",             "EList",
    "log2",        "limma",             "EList")

  method_params <- tribble(
    ~dge_method,    ~robust_fit,  ~robust_ebayes,  ~trend_ebayes,
    "voom",         FALSE,        FALSE,           FALSE,
    "edgeR-qlf",    TRUE,         FALSE,           FALSE,
    "limma-trend",  FALSE,        FALSE,           TRUE,
    "limma",        FALSE,        FALSE,           FALSE)

  info <- left_join(assay_methods, method_params, by = "dge_method")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}
