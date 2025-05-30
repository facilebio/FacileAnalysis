#' Peform a differential expression analysis.
#'
#' Use [flm_def()] to define the design matrix and contrast to test and
#' pass the `FacileLinearModelDefinition` object returned from that to `fdge()`
#' to run the desired differential testing framework (dictated by the `method`
#' parameter) over the data. `flm_def` accepts a
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
#'    mean, uses  [limma::arrayWeights()] when `with_sample_weights = TRUE`.
#' * `"limma"`: Straightup limma, this expects log2-normal like data, with
#'    (largely) no trend of the variance to the mean worth modeling. Uses
#'    [limma::arrayWeights()] when `with_sample_weights = TRUE`
#'
#' @section Feature Filtering Strategy:
#' You will almost always want to filter out lowly abundant features before
#' performing differential expression analysis. You can either do this by
#' explicitly requesting which features to test via the `features` parameter,
#' or by setting `filter = "default`.
#'
#' When `filter == "default"`, the filtering strategy is largely based on the
#' logic found in [edgeR::filterByExpr()].
#'
#' When `fdge` analysis is performed on count data, the filtering is precisely
#' executed using this function, using `design(x)` as the design parameter to
#' `filterByExpr`. You can modify the filtering behavior by passing some named
#' parameters found in the [edgeR::filterByExpr()] function down to it via
#' `fdge`'s `...` parameter, with the exception of the `design` paremeter,
#' because this is already defined. Use: 
#' * `filter_min_count` for the `edgeR::filterByExpr(min.count)` parameter
#' * `filter_min_total_count` for `edgeR::filterByExpr(min.total.count)`
#' * `filter_large_n` for `edgeR::filterByExpr(large.n)`; and 
#' * `filter_min_prop` for `edgeR::filterByExpr(min.prop)`.
#'
#' There are times when you want to tweak this behavior in ways that aren't
#' exactly supported by `filterByExpr`. You can pass in a "feature descriptor"
#' (a character vector of feature ids, or a data.frame with a "feature_id"
#' column) into the following parameters:
#'
#' * `filter_universe`: The features enumerated in this parameter will restrict
#'   the universe of features that can potentially be included in the downstream
#'   analysis. The `filterByExpr()` logic will happen downstream of this
#'   universe. The default value is `NULL`, which specifies the universe of
#'   features to be all of the ones measured using this assay.
#' * `filter_require`: The `filterByExpr` logic happens first on the universe
#'   of features as parameterized. All features enumerated here will be forcibly
#'   included in the analysis, irrespective of whether they would have passed
#'   the perscribed filter criteria or not. The defalut value for this argument
#'   is `NULL`, which means no genes are forcibly included in the analysis when
#'   they do not pass muster given the filtering criteria.
#'
#' @export
#' @param x a data source
#' @param assay_name the name of the assay that holds the measurements for test.
#'   Defaults to `default_assay(x)`.
#' @param method The differential testing framework to use over the data. The
#'   valid choices are defined by the type of assay `assay_name`is. Refer to the
#'   *Differential Expression Testing Methods* section for more details.
#' @param features Explicitly request the features to test. If this is provided,
#'   then the `filter` criteria (to remove low abundance features, for instance)
#'   is skipped. By default this is `NULL`
#' @param filter,with_sample_weights Passed into [biocbox()] to determine which
#'   features (genes) are removed from the dataset for testing, as well as if
#'   to use [limma::arrayWeights()] or [limma::voomWithQualityWeights()]
#'   (where appropriate) when testing (default is not to).
#' @param weights a `sample_id,feature_id,weight` table of observation weights
#'   to use when `method == "limma"`.
#' @param treat_lfc If this is numeric, this activates limma's "treat"
#'   functionality and tests for differential expression against this specified 
#'   log2 fold change threshold.
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#' samples <- efds |>
#'   FacileData::filter_samples(indication == "BLCA") |>
#'   FacileData::with_sample_covariates() |> 
#'   dplyr::mutate(something = sample(c("a", "b"), dplyr::n(), replace = TRUE))
#' mdef <- flm_def(samples, covariate = "sample_type",
#'                 numer = "tumor", denom = "normal",
#'                 batch = "sex",
#'                 metadata = list(label = "test flm"))
#' dge <- fdge(mdef, method = "voom", metadata = list(label = "test dge"))
#' if (interactive()) {
#'   viz(dge)
#'   viz(dge, "146909")
#'   shine(dge)
#' }
#' dge.stats <- tidy(dge)
#' dge.sig <- signature(dge)
#'
#' stage.anova <- samples |>
#'   flm_def(covariate = "stage", batch = "sex") |>
#'   fdge(method = "voom")
#' anova.sig <- signature(stage.anova)
fdge <- function(x, assay_name = NULL, method = NULL, features = NULL,
                 filter = "default", with_sample_weights = FALSE,
                 treat_lfc = NULL, ...,
                 metadata = list(),
                 verbose = FALSE) {
  UseMethod("fdge", x)
}

#' @noRd
#' @export
fdge.FacileFailedModelDefinition <- function(x, assay_name = NULL,
                                             method = NULL, features = NULL,
                                             filter = "default",
                                             with_sample_weights = FALSE, ...,
                                             metadata = list(),
                                             verbose = FALSE) {
  stop("There are erros in this model:\n", paste(x$errors, collapse = "\n"))
}


#' @export
#' @rdname fdge
fdge.FacileAnovaModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                            features = NULL, filter = "default",
                                            with_sample_weights = FALSE, ...,
                                            metadata = list(),
                                            verbose = FALSE) {
  res <- NextMethod(coef = x[["coef"]])
  # rename .intercept. to mean.level
  covariate <- param(x, "covariate")
  vals <- samples(x)[[covariate]]
  if (is.factor(vals)) {
    intercept <- levels(vals)[1L]
  } else {
    intercept <- sort(unique(vals))[1L]
  }
  intercept <- paste0("mean.", intercept)
  data.table::setnames(res$result, ".intercept.", intercept)
  res
}

#' @export
#' @rdname fdge
fdge.FacileTtestDGEModelDefinition <- function(x, assay_name = NULL,
                                               method = NULL, features = NULL,
                                               filter = "default",
                                               with_sample_weights = FALSE,
                                               treat_lfc = NULL,
                                               ...,
                                               metadata = list(),
                                               verbose = FALSE) {
  res <- NextMethod(contrast = x[["contrast"]])
  res
}


#' @export
#' @rdname fdge
#' @param ... passed down into inner methods, such as `biocbox` to tweak
#'   filtering criteria, for instance
fdge.FacileLinearModelDefinition <- function(x, assay_name = NULL,
                                             method = NULL, features = NULL,
                                             filter = "default",
                                             with_sample_weights = FALSE,
                                             treat_lfc = NULL,
                                             weights = NULL, with_box = FALSE,
                                             ...,
                                             biocbox = NULL,
                                             trend.eBayes = FALSE,
                                             robust.eBayes = FALSE,
                                             metadata = list(),
                                             verbose = FALSE) {
  messages <- character()
  warnings <- character()
  errors <- character()
  .fds <- assert_facile_data_store(fds(x))

  if (is.null(assay_name)) assay_name <- default_assay(.fds)
  
  ainfo <- try(assay_info(.fds, assay_name), silent = TRUE)
  bad_assay <- is(ainfo, "try-error")

  if (bad_assay) {
    msg <- glue("assay_name `{assay_name}` not present in FacileDataStore")
    errors <- c(errors, msg)
    assay_type <- "ERROR"
  } else {
    assay_type <- ainfo$assay_type
  }
  
  if (bad_assay) {
    dge_methods <- filter(fdge_methods(), FALSE)
    dropped <- distinct(samples(x), dataset, sample_id, .keep_all = TRUE)
    assay_samples <- filter(dropped, FALSE)
  } else {
    dge_methods <- fdge_methods(assay_type)
    assay_samples <- filter_by_assay_support(samples(x), assay_name)
    dropped <- samples(assay_samples, dropped = TRUE)
  }
    
  if (nrow(dge_methods) == 0L) {
    msg <- glue("Differential expression method `{method}` for assay_type ",
                "`{assay_type}` not found")
    errors <- c(errors, msg)
  }
  
  ndropped <- nrow(dropped)
  if (ndropped > 0) {
    xo <- x
    msg <- paste(
      ndropped, "samples have no", assay_name, "data. These samples will ",
      "be removed for downstream analysis.")
    ftrace(msg)
    
    # Re-estimate the linear model
    # NOTE: This should be moved to biocbox.FacileLinearModelDefinition
    if (length(errors) == 0L) {
      msg <- paste(
        "re-estimating linear model from reduced sample space",
        "from", nrow(samples(xo)), "->", nrow(assay_samples))
      ftrace(msg)
      x <- redo(x, samples = assay_samples)
      errors <- c(errors, errors(x))
      warnings <- c(warnings, warnings(x))
      messages <- c(msg, messages(x))
    }
  }
  
  if (length(errors) == 0L) {
    if (!is.null(biocbox)) {
      # We are doing this so we don't have to recalculate things like
      # edgeR::estimateDisp, or limma::voom
      ftrace("reusing a cached biocbox")
      messages <- c(messages, "using a cached biobox needs more testing")
      if (is.null(method)) {
        stop("Using a cached biocbox requires you to explicitly set `method` ",
             "parameter")
      }
      bb <- biocbox
      fbits <- attr(bb, "facile")

      if (method != fbits[["params"]][["method"]]) {
        msg <- paste0(
          "The previously used `method` ('",
          fbits[["params"]][["method"]],
          "'), does not match the one currently being run ('", method, "')")
        warnings <- c(warnings, msg)
        fbits[["params"]][["method"]] <- method
      }
      
      if (!is.null(weights)) {
        bb <- .add_observation_weights(bb, weights)
      }
      
      messages <- c(messages, fbits[["messages"]])
      warnings <- c(warnings, fbits[["warnings"]])
      errors <- c(errors, fbits[["errors"]])
    } else {
      ftrace("... retrieving expression data")
      
      bb <- biocbox(x, assay_name, method, features, filter,
                    with_sample_weights = with_sample_weights,
                    weights = weights, ...)
      fbits <- attr(bb, "facile")
      messages <- c(messages, fbits[["messages"]])
      warnings <- c(warnings, fbits[["warnings"]])
      errors <- c(errors, fbits[["errors"]])
    }

    method <- fbits[["params"]][["method"]]

    des <- design(x)
    if (!setequal(rownames(des), colnames(bb))) {
      errors <- c(errors, "rownames of design(x) does not match biocbox")
      # return(out)
      stop(errors)
    }

    if (!isTRUE(all.equal(rownames(des), colnames(bb)))) {
      bb <- bb[,rownames(des)]
    }

    if (!is.null(treat_lfc)) {
      if (!test_number(treat_lfc, lower = 0) || !is_ttest(x)) {
        warnings <- c(
          warnings,
          "Illegal parameter passed to `treat_lfc`. It is being ignored")
        treat_lfc <- 0
      }
    } else {
      treat_lfc <- 0
    }
    use.treat <- treat_lfc > 0

    if (verbose) {
      message("... running differential expression analysis")
    }

    # TODO: Add more params to send to calculateIndividualLogFC based on
    # dge method, for instance when for:
    # * method == "limma-trend", send trend.eBayes = TRUE
    if (method == "limma-trend") {
      trend.eBayes <- TRUE
    }

    block <- param(x, "block")
    dup.corr <- NULL
    if (!is.null(block)) {
      # this needs to be a EList for now
      assert_class(bb, "EList")
      block <- bb[["targets"]][[block]]
      assert_categorical(block, len = ncol(bb), any.missing = FALSE)
      dup.corr <- bb$block.corr
    }

    if (is_ttest(x)) {
      testme <- x[["contrast"]]
      clazz <- "FacileTtestAnalysisResult"
    } else {
      testme <- x[["coef"]]
      clazz <- "FacileAnovaAnalysisResult"
    }
    
    result <- sparrow::calculateIndividualLogFC(
      bb, des, contrast = testme,
      treat.lfc = treat_lfc,
      trend.eBayes = trend.eBayes,
      robust.eBayes = robust.eBayes,
      # with.fit = TRUE,
      block = block, correlation = dup.corr)

    # sparrow::calculateIndividualLogFC returns the stats table ordered by
    # featureId, let's put the features back in the order they are in y
    rownames(result) <- result[["feature_id"]]
    result <- result[rownames(bb),]

    axe.cols <- c("featureId", "x.idx")
    if (is(bb, "DGEList")) {
      axe.cols <- c(axe.cols, "t")
    }
    for (col in axe.cols) {
      if (col %in% names(result)) result[[col]] <- NULL
    }
    result <- as_tibble(result)
  } else {
    result <- NULL
  }

  # Hack here to support call from a bioc-container that we haven't turned
  # into a facile data store yet
  if (!is.null(result) && !"feature_type" %in% colnames(result)) {
    feature_type <- FacileData::infer_feature_type(
      result[["feature_id"]], summarize = TRUE)
    result[["feature_type"]] <- feature_type[["feature_type"]]
    result <- select(result, feature_type, everything())
  }

  # We'll stick on the lib.size and norm.factors on these samples so that
  # downstream queries over this result will use the same normalization factors
  # when retrieving data as well.
  if (is(bb, "DGEList") || is(bb, "EList")) {
    samples. <- samples(x)
    sinfo <- if (is(bb, "DGEList")) bb[["samples"]] else bb[["targets"]]
    for (cname in c("lib.size", "norm.factors")) {
      if (cname %in% colnames(sinfo) && is.numeric(sinfo[[cname]])) {
        samples.[[cname]] <- sinfo[[cname]]
      }
    }
    x[["covariates"]] <- samples. # `samples(x) <- samples.` would be nice!
  }

  out <- list(
    biocbox = bb,
    result = result,
    params = list(
      assay_name = assay_name,
      method = method,
      model_def = x,
      treat_lfc = if (use.treat) treat_lfc else NULL,
      with_sample_weights = with_sample_weights,
      weights = weights),
    # Standard FacileAnalysisResult things
    fds = .fds,
    metadata = metadata,
    messages = messages,
    warnings = warnings,
    errors = errors)

  if (with_box) {
    out[["biocbox"]] <- bb
  }
  class(out) <- c(clazz, "FacileDgeAnalysisResult", "FacileAnalysisResult")
  out
}

# Methods and Accessors ========================================================

#' @noRd
#' @export
result.FacileDgeAnalysisResult <- function(x, name = "result", ...) {
  tidy(x, name, ...)
}

#' @param features you can provide a set of features to reduce the analysis
#'   to. This will readjust the `padj` values to the features used.
#' @param padjust The method used to readjust pvalues. Passed to `method`
#'   parameter in [stats::p.adjust()].
#' @noRd
#' @export
tidy.FacileDgeAnalysisResult <- function(x, name = "result", features = NULL,
                                         padjust = "BH", ...) {
  out <- x[["result"]]
  if (!is.null(features)) {
    features <- extract_feature_id(features, as_tibble = TRUE)
    out <- left_join(features, out, by = "feature_id")
    if (nrow(out) == 0L) {
      stop("None of the requested features were found in the resut.")
    }
    missed <- anti_join(features, out, by = "feature_id")
    nmiss <- nrow(missed)
    if (nmiss) {
      warning(length(nmiss), " features not found in result")
    }
    out[["padj"]] <- stats::p.adjust(out[["pval"]], padjust)
  }
  out
}

#' @noRd
#' @export
features.FacileDgeAnalysisResult <- function(x, ...) {
  stat.table <- tidy(x)
  take <- c("feature_id", "feature_type", "symbol", "assay", "assay_type")
  take <- intersect(take, colnames(stat.table))
  select(stat.table, {{take}}, everything())
}

#' @noRd
#' @export
model.FacileDgeAnalysisResult <- function(x, ...) {
  param(x, "model_def")
}

#' @noRd
#' @export
design.FacileDgeAnalysisResult <- function(x, ...) {
  design(model(x, ...), ...)
}

#' @noRd
#' @export
contrast.FacileDgeAnalysisResult <- function(x, ...) {
  contrast(param(x, "model_def"), ...)
}

#' @noRd
#' @export
name.FacileTtestAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  batch <- param(m, "batch")
  contrast <- gsub(" ", "", m[["contrast_string"]])
  contrast <- gsub("/", ".divby.", contrast, fixed = TRUE)
  out <- sprintf("%s_%s", covname, contrast)
  if (length(batch)) {
    out <- sprintf("%s_controlfor_%s", out, paste(batch, collapse = ","))
  }
  out
}

#' @noRd
#' @export
label.FacileTtestAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  batch <- param(m, "batch")
  contrast <- m[["contrast_string"]]
  out <- sprintf("%s: %s", covname, contrast)
  if (length(batch)) {
    out <- sprintf("%s (control for %s)", out, paste(batch, collapse = ","))
  }
  out
}

#' @noRd
#' @export
name.FacileAnovaAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  clevels <- unique(samples(m)[[covname]])
  batch <- param(m, "batch")
  out <- sprintf("%s_%s", covname, paste(clevels, collapse = ","))
  if (length(batch)) {
    out <- sprintf("%s_controlfor_%s", out, paste(batch, collapse = ","))
  }
  out
}

#' @noRd
#' @export
label.FacileAnovaAnalysisResult <- function(x, ...) {
  m <- model(x)
  covname <- param(m, "covariate")
  clevels <- unique(samples(m)[[covname]])
  batch <- param(m, "batch")
  out <- sprintf("%s [%s]", covname, paste(clevels, collapse = ","))
  if (length(batch)) {
    out <- sprintf("%s (control for %s)", out, paste(batch, collapse = ","))
  }
  out
}

# Ttest Ranks and Signatures ===================================================

#' @export
#' @noRd
ranks.FacileTtestAnalysisResult <- function(x, signed = TRUE, rank_by = "logFC",
                                            ...) {
  ranks. <- tidy(x, ...)
  rank_by <- assert_choice(rank_by, colnames(ranks.))
  assert_numeric(ranks.[[rank_by]])
  ofn <- if (rank_by %in% c("pval", "padj")) identity else dplyr::desc
  if (signed) {
    ranks. <- arrange(ranks., ofn(ranks.[[rank_by]]))
  } else {
    ranks. <- arrange(ranks., ofn(abs(ranks.[[rank_by]])))
  }

  out <- list(
    result = ranks.,
    params = list(signed = signed),
    ranking_columns = "logFC",
    ranking_order = "descending",
    fds = fds(x))

  # FacileTtestFeatureRanksSigned
  clazz <- "Facile%sFeatureRanks%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c("Ttest", "Ttest",  "", ""), c(s, "", s, ""))
  class(out) <- c(classes, "FacileFeatureRanks")
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
signature.FacileTtestFeatureRanks <- function(x, min_logFC = x[["treat_lfc"]],
                                              max_padj = 0.10,
                                              ntop = 20,
                                              name = NULL,
                                              collection_name = NULL,
                                              ...) {
  if (is.null(name)) name <- "Ttest signature"
  name. <- assert_string(name)

  if (is.null(collection_name)) collection_name <- class(x)[1L]
  assert_string(collection_name)

  if (is.null(min_logFC)) min_logFC <- 0
  min_logFC <- abs(min_logFC)

  res <- tidy(x) |>
    mutate(direction = ifelse(logFC > 0, "up", "down"))

  if (signed(x)) {
    up <- res |>
      filter(padj <= max_padj, logFC > min_logFC) |>
      head(ntop)
    down <- res |>
      filter(padj <= max_padj, logFC < -min_logFC) |>
      tail(ntop)
    sig <- bind_rows(up, down)
  } else {
    sig <- res |>
      filter(padj <= max_padj) |>
      head(ntop)
  }

  sig <- sig |>
    mutate(collection = collection_name, name = paste(name., direction)) |>
    select(collection, name, feature_id, symbol, direction,
           logFC, pval, padj, everything())

  out <- list(
    result = sig,
    params = list(max_padj = max_padj, min_logFC = min_logFC, ntop = ntop))
  class(out) <- sub("Ranks", "Signature", class(x))
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
  ranks. <- arrange(result(x, ...), desc(.data[["F"]]))

  out <- list(
    result = ranks.,
    params = list(signed = signed),
    ranking_columns = "F",
    ranking_order = "descending")
  # FacileAnovaFeatureRanksSigned
  clazz <- "Facile%sFeatureRanks%s"
  s <- if (signed) "Signed" else "Unsigned"
  classes <- sprintf(clazz, c("Anova", "Anova",  "", ""), c(s, "", s, ""))
  class(out) <- c(classes, "FacileFeatureRanks")
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
signature.FacileAnovaFeatureRanks <- function(x, max_padj = 0.10,
                                              ntop = 20, name = NULL,
                                              collection_name = NULL,
                                              ...) {
  if (is.null(name)) name <- "Anova signature"
  name. <- assert_string(name)

  if (is.null(collection_name)) collection_name <- class(x)[1L]
  assert_string(collection_name)

  sig <- result(x) |>
    filter(padj <= max_padj) |>
    mutate(collection = collection_name, name = name.,
           direction = "udefined") |>
    head(ntop) |>
    select(collection, name, feature_id, symbol, direction, F, pval, padj,
           everything())

  out <- list(
    result = sig,
    params = list(ntop = ntop, max_padj = max_padj))
  class(out) <- sub("Ranks$", "Signature", class(x))
  out
}

# Facile API ===================================================================

#' @noRd
#' @export
samples.FacileDgeAnalysisResult <- function(x, ...) {
  samples(model(x), ...)
}

#' @section FacileDgeAnalysisResult:
#' Given a FacileDgeAnalysisResult, we can re-materialize the Bioconductor assay
#' container used within the differential testing pipeline used from [fdge()].
#' Currently we have limited our analysis framework to either work over DGEList
#' (edgeR) or EList (limma) containers.
#'
#' @export
#' @rdname biocbox
biocbox.FacileDgeAnalysisResult <- function(x, cached = TRUE, ...) {
  if (cached) {
    out <- x[["biocbox"]]
  } else {
    # I originally had this method signature as (x, ...) and then explicitly
    # passed down assay_name = param(x, "assay_name"), but if we don't catch
    # this arguments in the function and the user calls the function with them,
    # we get an error of duplicate parameters in the function call when we
    # delegate down to biocbox.FacileLinearModelDefinition
    assay_name <- assert_string(param(x, "assay_name"))
    method <- assert_string(param(x, "method"))

    # meta.cols <- c("lib.size", "norm.factors", "sizeFactor")
    # meta.cols <- intersect(meta.cols, )
    out <- biocbox(model(x), assay_name = assay_name, method = method,
                   features = features(x),
                   with_sample_weights = param(x, "with_sample_weights"),
                   weights = param(x, "weights"), ...)
  }
  out
}

#' @noRd
#' @export
print.FacileDgeAnalysisResult <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
format.FacileDgeAnalysisResult <- function(x, ...) {
  test_type <- if (is_ttest(x)) "ttest" else "ANOVA"
  mdef <- model(x)
  des <- design(mdef)
  formula <- mdef[["design_formula"]]

  if (test_type == "ttest") {
    test <- mdef[["contrast_string"]]
    des.cols <- names(mdef$contrast[mdef$contrast != 0])
  } else {
    test <- sprintf("%s (%s)", param(mdef, "covariate"),
                    paste(colnames(des)[mdef[["coef"]]], collapse = "|"))
    des.cols <- c(1, mdef$coef)
  }

  # nsamples <- sum(colSums(des[, des.cols, drop = FALSE]) != 0)
  nsamples <- sum(rowSums(des[,des.cols,drop = FALSE]) > 0)

  res <- result(x)
  ntested <- nrow(res)
  nsig <- sum(!is.na(res[["padj"]]) & res[["padj"]] < 0.10)
  nsig2 <- sum(!is.na(res[["padj"]]) & res[["padj"]] < 0.10 & abs(res[["logFC"]]) > 1)

  out <- paste(
    "=======================================================================\n",
    sprintf("FacileDgeAnalysisResult (%s)\n", test_type),
    "-----------------------------------------------------------------------\n",
    glue(
      "Significant Results (FDR < 0.1): ",
      "{nsig} / {ntested} [{nsig2} > 2fold]", .trim = FALSE), "\n",
      "Formula: ", formula, "\n",
      "Tested: ", test, "\n",
      "Number of samples: ", nsamples, "\n",
    "=======================================================================\n",
    sep = "")
  out
}

# Helpers ======================================================================


#' Defines a voomLmFit method that can accept dots (and toss them)
#'
#' This can't work in the current framework because these methods need to
#' return an EList (biocbox), not a fit object
#'
#' importFrom edgeR voomLmFit
#' @noRd
.voomLmFit <- function(counts, design = NULL, block = NULL,
                       prior.weights = NULL, sample.weights = FALSE, ...) {
  args <- list(...)
  vm.args <- formals(edgeR::voomLmFit)
  take.args <- intersect(names(vm.args), names(args))

  call.args <- list(counts = counts, design = design, block = block,
                    prior.weights = prior.weights,
                    sample.weights = sample.weights)
  if (length(take.args)) call.args <- c(call.args, args[take.args])

  do.call(edgeR::voomLmFit, call.args)
}


#' Defines a voom method that can accept dots (and toss them)
#'
#' @noRd
#' @importFrom limma voom duplicateCorrelation
.voom <- function(counts, design = NULL, block = NULL, save.plot = TRUE, ...) {
  args <- list(...)
  vm.args <- formals(limma::voom)
  take.args <- intersect(names(vm.args), names(args))

  call.args <- list(counts = counts, design = design, save.plot = save.plot)
  if (length(take.args)) call.args <- c(call.args, args[take.args])

  vm <- do.call(limma::voom, call.args)

  # Follow the two-step dupcor approach outlined in Section 18.1.9 of the
  # limmaUserGuide
  if (!is.null(block)) {
    dcor <- duplicateCorrelation(vm, vm$design, block = vm$targets[[block]])

    call.args[["correlation"]] <- dcor$consensus.correlation
    call.args[["block"]] <- vm$targets[[block]]

    vm <- do.call(limma::voom, call.args)
    dcor <- duplicateCorrelation(vm, vm$design, block = vm$targets[[block]])
    vm$block.corr <- dcor$consensus.correlation
  }

  vm
}

#' This is necessary because limma::voomWithQualityWeights calls voom with
#' voom(...), but it doesn't accept ...
#' @noRd
#' @importFrom limma voomWithQualityWeights duplicateCorrelation
.voomWithQualityWeights <- function(counts, design = NULL, block = block, 
                                    save.plot = TRUE, ...) {
  args <- list(...)
  vm.args <- formals(limma::voom)
  vmw.args <- formals(limma::voomWithQualityWeights)
  all.formals <- c(vm.args, vmw.args)
  take.args <- intersect(names(args), names(all.formals))

  call.args <- list(counts = counts, design = design, save.plot = save.plot)
  if (length(take.args)) call.args <- c(call.args, args[take.args])

  vm <- do.call(limma::voomWithQualityWeights, call.args)

  # Follow the two-step dupcor approach outlined in Section 18.1.9 of the
  # limmaUserGuide
  if (!is.null(block)) {
    dcor <- duplicateCorrelation(vm, vm$design, block = vm$targets[[block]])

    call.args[["correlation"]] <- dcor$consensus
    call.args[["block"]] <- vm$targets[[block]]

    vm <- do.call(limma::voom, call.args)
    dcor <- duplicateCorrelation(vm, vm$design, block = vm$targets[[block]])
    vm$block.corr <- dcor$consensus
  }

  vm
}

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
fdge_methods <- function(assay_type = NULL,
                         on_missing = c("error", "warning")) {
  on_missing <- match.arg(on_missing)
  # assay_type values : rnaseq, umi, affymrna, affymirna, log2

  # This is a table of assay_type : dge_method possibilites. The first row
  # for each assay_type is the default analysis method
  assay_methods <- read.csv(
    system.file("extdata", "analysis-params", "fdge", "assay-methods.csv",
                package = "FacileAnalysis"),
    strip.white = TRUE)

  method_params <- read.csv(
    system.file("extdata", "analysis-params", "fdge", "method-params.csv",
                package = "FacileAnalysis"),
    strip.white = TRUE)

  info <- left_join(assay_methods, method_params, by = "dge_method")

  if (!is.null(assay_type)) {
    if (on_missing == "error") {
      assert_choice(assay_type, info[["assay_type"]])
    }
    info <- filter(info, .data$assay_type == .env$assay_type)
  }

  info
}
