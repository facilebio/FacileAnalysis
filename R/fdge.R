#' Peform a differential expression analysis.
#'
#' Use [fdge_model_def()] to define the design matrix and contrast to test and
#' pass the `FacileDgeModelDefinition` object returned from that to `fdge()`
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
#' @importFrom multiGSEA calculateIndividualLogFC logFC multiGSEA
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
#'                        numer = "tumor", denom = "normal", fixed = "sex")
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
#' @param ... passed down into inner methods, such as `biocbox` to tweak
#'   filtering criteria, for instance
fdge.FacileDgeModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                          filter = "default",
                                          with_sample_weights = FALSE,
                                          treat_lfc = NULL,
                                          ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  .fds <- assert_facile_data_store(fds(x))

  if (is.ttest(x)) {
    testme <- x[["contrast"]]
    clazz <- "FacileTtestAnalysisResult"
  } else {
    testme <- x[["coef"]]
    clazz <- "FacileAnovaAnalysisResult"
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

    if (!is.null(treat_lfc)) {
      if (!test_number(treat_lfc, lower = 0) || !is.ttest(x)) {
        warnings <- c(
          warnings,
          "Illegal parameter passed to `treat_lfc`. It is being ignored")
        treat_lfc <- 0
      }
    } else {
      treat_lfc <- 0
    }
    use.treat <- treat_lfc > 0
    y <- result(bb)
    des <- design(bb)
    result <- calculateIndividualLogFC(y, des, contrast = testme,
                                       treat.lfc = treat_lfc)
    # multiGSEA::calculateIndividualLogFC returns the stats table ordered by
    # featureId, let's put the features back in the order they are in y
    rownames(result) <- result[["featureId"]]
    result <- result[rownames(y),]

    axe.cols <- c("featureId", "x.idx")
    if (is(y, "DGEList")) {
      axe.cols <- c(axe.cols, "t")
    }
    for (col in axe.cols) {
      if (col %in% names(result)) result[[col]] <- NULL
    }
    result <- as.tbl(result)
  } else {
    result <- NULL
  }

  out <- list(
    result = result,
    params = list(
      assay_name = assay_name,
      method = bb[["dge_method"]],
      model_def = x,
      treat_lfc = if (use.treat) treat_lfc else NULL,
      with_sample_weights = with_sample_weights),
    # Standard FacileAnalysisResult things
    fds = .fds,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(out) <- c(clazz, "FacileDgeAnalysisResult", "FacileAnalysisResult")
  out
}

# Methods and Accessors ========================================================

#' @noRd
#' @export
result.FacileDgeAnalysisResult <- function(x, name, ...) {
  out <- x[["result"]]
  if (!"feature_type" %in% colnames(out)) {
    feature_type <- guess_feature_type(out[["feature_id"]], summarize = TRUE)
    out[["feature_type"]] <- feature_type[["feature_type"]]
    out <- select(feature_type, everything())
  }
  out
}

#' @noRd
#' @export
model.FacileDgeAnalysisResult <- function(x, ...) {
  x[["params"]][["model_def"]]
}

#' @noRd
#' @export
design.FacileDgeAnalysisResult <- function(x, ...) {
  design(model(x, ...), ...)
}

#' @noRd
#' @export
name.FacileTtestAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  fixed <- param(m, "fixed")
  contrast <- gsub(" ", "", m[["contrast_string"]])
  contrast <- gsub("/", ".divby.", contrast, fixed = TRUE)
  out <- sprintf("%s_%s", covname, contrast)
  if (length(fixed)) {
    out <- sprintf("%s_controlfor_%s", out, paste(fixed, collapse = ","))
  }
  out
}

#' @noRd
#' @export
label.FacileTtestAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  fixed <- param(m, "fixed")
  contrast <- m[["contrast_string"]]
  out <- sprintf("%s: %s", covname, contrast)
  if (length(fixed)) {
    out <- sprintf("%s (control for %s)", out, paste(fixed, collapse = ","))
  }
  out
}

#' @noRd
#' @export
name.FacileAnovaAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  clevels <- unique(samples(m)[[covname]])
  fixed <- param(m, "fixed")
  out <- sprintf("%s_%s", covname, paste(clevels, collapse = ","))
  if (length(fixed)) {
    out <- sprintf("%s_controlfor_%s", out, paste(fixed, collapse = ","))
  }
  out
}

#' @noRd
#' @export
label.FacileAnovaAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  clevels <- unique(samples(m)[[covname]])
  fixed <- param(m, "fixed")
  out <- sprintf("%s [%s]", covname, paste(clevels, collapse = ","))
  if (length(fixed)) {
    out <- sprintf("%s (control for %s)", out, paste(fixed, collapse = ","))
  }
  out
}

# Ttest Ranks and Signatures ===================================================

#' @export
#' @noRd
ranks.FacileTtestAnalysisResult <- function(x, signed = TRUE, ...) {
  ranks. <- result(x, ...)
  if (signed) {
    ranks. <- arrange(ranks., desc(logFC))
  } else {
    ranks. <- arrange(ranks., pval)
  }

  out <- list(
    result = ranks.,
    ranking_columns = "logFC",
    ranking_order = "descending")
  # FacileTtestFeatureRanksSigned
  clazz <- "Facile%sFeatureRanks%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c("Ttest", "Ttest",  "", ""), c(s, "", s, ""))
  class(out) <- c(clazz,
                  "FacileFeatureRanks",
                  "FacileAnalysisResult")
  out
}

#' @export
#' @noRd
signature.FacileTtestAnalysisResult <- function(x, min_logFC = x[["treat_lfc"]],
                                                max_padj = 0.10, ntop = 20,
                                                name = NULL,
                                                collection_name = NULL,
                                                ...) {
  signature(ranks(x, ...), min_logFC = min_logFC, max_padj = max_padj,
            ntop = ntop, name = name, collection_name = collection_name, ...)
}

#' @export
#' @noRd
signature.FacileTtestFeatureRankings <- function(x, min_logFC = x[["treat_lfc"]],
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

# ANOVA Ranks and Signatures ===================================================

#' @noRd
#' @export
ranks.FacileAnovaAnalysisResult <- function(x, signed = FALSE, ...) {
  if (signed) {
    warning("ANOVA results can only provide unsigned ranks ",
            "based on p-value, which is the same as desc(Fstatistic)",
            immediate. = TRUE)
    signed <- FALSE
  }
  ranks. <- arrange(result(x, ...), desc(F))

  out <- list(
    result = ranks.,
    ranking_columns = "F",
    ranking_order = "descending")
  # FacileAnovaFeatureRanksSigned
  clazz <- "Facile%sFeatureRanks%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c("Anova", "Anova",  "", ""), c(s, "", s, ""))
  class(out) <- c(clazz,
                  "FacileFeatureRanks",
                  "FacileAnalysisResult")
  out
}

#' @export
#' @noRd
signature.FacileAnovaAnalysisResult <- function(x, max_padj = 0.10, ntop = 20,
                                           name = NULL, collection_name = NULL,
                                           ...) {
  signature(ranks(x, ...), max_padj = max_padj, ntop = ntop, name = name,
            collectoin_name = collection_name, ...)
}

#' @export
#' @noRd
signature.FacileAnovaFeatureRankings <- function(x, max_padj = 0.10,
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
samples.FacileDgeAnalysisResult <- function(x, ...) {
  samples(model(x))
}

#' @section FacileDgeAnalysisResult:
#' Given a FacileDgeAnalysisResult, we can re-materialize the Bioconductor assay
#' container used within the differential testing pipeline used from [fdge()].
#' Currently we have limited our analysis framework to either work over DGEList
#' (edgeR) or EList (limma) containers.
#'
#' @export
#' @rdname biocbox
biocbox.FacileDgeAnalysisResult <- function(x, ...) {
  res <- biocbox(model(x),
                 x[["params"]][["assay_name"]],
                 x[["params"]][["method"]],
                 filter = result(x)[["feature_id"]], ...)
  res
}

#' @noRd
#' @export
print.FacileDgeAnalysisResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.FacileDgeAnalysisResult <- function(x, ...) {
  test_type <- if (is.ttest(x)) "ttest" else "ANOVA"
  mdef <- model(x)
  des <- design(mdef)
  formula <- mdef[["design_formula"]]

  if (test_type == "ttest") {
    test <- mdef[["contrast_string"]]
    des.cols <- names(mdef$contrast[mdef$contrast != 0])
  } else {
    test <- sprintf("%s (%s)", mdef[["covariate"]],
                    paste(colnames(des)[mdef[["coef"]]], collapse = "|"))
    des.cols <- c(1, mdef$coef)
  }

  nsamples <- sum(colSums(des[, des.cols, drop = FALSE] != 0))

  res <- result(x)
  ntested <- nrow(res)
  nsig <- sum(!is.na(res[["padj"]]) & res[["padj"]] < 0.10)

  out <- paste(
    "===========================================================\n",
    sprintf("FacileDgeAnalysisResult (%s)\n", test_type),
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
    "cpm",         "limma-trend",       "EList",
    "isoseq",      "voom",              "EList",
    "isoseq",      "limma-trend",       "EList",
    "affymrna",    "limma",             "EList",
    "affymirna",   "limma",             "EList",
    "lognorm",     "limma",             "EList")

  method_params <- tribble(
    ~dge_method,    ~robust_fit,  ~robust_ebayes,  ~trend_ebayes, ~can_sample_weight,
    "voom",         FALSE,        FALSE,           FALSE,          TRUE,
    "edgeR-qlf",    TRUE,         FALSE,           FALSE,          FALSE,
    "limma-trend",  FALSE,        FALSE,           TRUE,           TRUE,
    "limma",        FALSE,        FALSE,           FALSE,          TRUE)

  info <- left_join(assay_methods, method_params, by = "dge_method")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}
