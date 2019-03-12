#' Peforms a differential expression analysis and GSEA.
#'
#' @export
#' @importFrom multiGSEA GeneSetDb calculateIndividualLogFC logFC multiGSEA
#'
#' @param x a data source
#' @param numer character vector defining the covariate/groups that
#'   make up the numerator
#' @param denom character vector defining the covariate/groups that
#'   make up the denominator
#' @param fixed character vector defining the covariate/groups to
#'   use as fixed effects
#' @param covariates a data.frame of covariates (columns) for each
#'   sample (rows) in `x`
#' @examples
#' mdef <- FacileData::exampleFacileDataSet() %>%
#'  filter_samples(indication == "BLCA") %>%
#'  fdge_model_def(covariate = "sample_type",
#'                 numer = "tumor",
#'                 denom = "normal")
#' dge <- fdge(mdef, method = "voom")
#'
#' if (interactive()) {
#'   viz(dge)
#'   shine(dge)
#' }
fdge <- function(x, ...) {
  UseMethod("fdge", x)
}

#' @export
#' @rdname fdge
fdge.FacileAnovaModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                            filter = "default", ...) {
  res <- NextMethod(coef = x$coef)
  res
}

#' @export
fdge.FacileTtestDGEModelDefinition <- function(x, assay_name = NULL,
                                               method = NULL,
                                               filter = "default",
                                               treat_lfc = NULL, ...) {
  res <- NextMethod(contrast = x$contrast)
  res
}

#' @export
#' @rdname fdge
fdge.FacileDGEModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                          filter = "default", treat_lfc = NULL,
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
    y <- biocbox(x, assay_name, method, dge_methods, filter, ...)
    messages <- c(messages, attr(y, "messages"))
    warnings <- c(warnings, attr(y, "warnings"))
    errors <- c(errors, attr(y, "errors"))

    method <- y$dge_method

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

    result <- calculateIndividualLogFC(y, y$design, contrast = testme,
                                       use.treat = use.treat,
                                       feature.min.logFC = treat_lfc)
    axe.cols <- c("featureId", "x.idx")
    for (col in axe.cols) {
      if (col %in% names(result)) result[[col]] <- NULL
    }
    result <- as.tbl(result)
    if (test_type == "ttest") {
      result <- arrange(result, desc(logFC))
    } else {
      result <- arrange(result, pval)
    }
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
samples.FacileDGEResult <- function(x, ...) {
  samples(x[["model_def"]])
}

#' @section FacileDGEResult:
#' Given a FacileDGEResult, we can re-materialize the Bioconductor assay
#' container used within the differential testing pipeline used from [fdge()].
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
#' @method print FacilePCAResult
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

.dge_methods <- c("voom", "qlf", "trended", "limma")

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
    ~assay_type, ~dge_method, ~bioc_class,
    "rnaseq",    "voom",      "DGEList",
    "rnaseq",    "qlf",       "DGEList",
    "rnaseq",    "trended",   "DGEList",
    "umi",       "voom",      "DGEList",
    "umi",       "qlf",       "DGEList",
    "umi",       "trended",   "DGEList",
    "tpm",       "trended",   "EList",
    "affymrna",  "limma",     "EList",
    "affymirna", "limma",     "EList",
    "log2",      "limma",     "EList")

  method_params <- tribble(
    ~dge_method,  ~robust_fit,  ~robust_ebayes,  ~trend_ebayes,
    "voom",       FALSE,        FALSE,           FALSE,
    "qlf",        TRUE,         FALSE,           FALSE,
    "trended",    FALSE,        FALSE,           TRUE,
    "limma",      FALSE,        FALSE,           FALSE)

  info <- left_join(assay_methods, method_params, by = "dge_method")

  if (!is.null(assay_type)) {
    assert_choice(assay_type, info[["assay_type"]])
    info <- info[info[["assay_type"]] == assay_type,]
  }

  info
}

# Other Helpers ================================================================

#' @noRd
#' @export
tidy.FacileDGEResult <- function(x, result = "dge") {
  # result <- match.arg(result, c("dge", multiGSEA::resultNames(result$gsea)))
  result <- match.arg(result, c("dge", multiGSEA::resultNames(x$gsea)))
  if (result == "dge") {
    out <- x[["dge"]]
  } else {
    # out <- multiGSEA::result(x$gsea, result)
  }
  out
}
