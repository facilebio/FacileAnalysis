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
fdge <- function(x, ...) {
  UseMethod("fdge", x)
}

#' @export
#' @rdname fdge
fdge.FacileAnovaModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                            gsea = NULL, filter = "default",
                                            ...) {
  res <- NextMethod(coef = x$coef)
  res
}

#' @export
fdge.FacileTtestDGEModelDefinition <- function(x, assay_name = NULL,
                                               method = NULL, gsea = "cameraPR",
                                               filter = "default", ...) {
  res <- NextMethod(contrast = x$contrast)
  res
}

fdge.FacileDGEModelDefinition <- function(x, assay_name = NULL, method = NULL,
                                          gsea = NULL, filter = "default",
                                          ...) {
  messages <- character()
  warnings <- character()
  errors <- character()

  clazz <- switch(class(x)[1L],
                  FacileTtestDGEModelDefinition = "FacileTtestDGEResult",
                  FacileAnovaModelDefinition = "FacileAnovaDGEResult",
                  NULL)

  .fds <- assert_class(x$fds, "FacileDataStore")

  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
  }

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

    # Do DGE and GSEA
    testme <- if (is.null(x[["coef"]])) x[["contrast"]] else x[["coef"]]

    # dge <- calculateIndividualLogFC(y, y$design, contrast = testme)
    gdb <- GeneSetDb(list(dummy = list(dummy = head(rownames(y), 5))))
    mg <- multiGSEA(gdb, y, y$design, contrast = testme, methods = "logFC")
  } else {
    mg <- NULL
  }

  out <- list(
    result = mg,
    assay_name = assay_name,
    method = method,
    model_def = x,
    # Standard FacileAnalysisResult things
    fds = .fds,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(out) <- c(clazz, "FacileDGEResult", "FacileAnalysisResult")
  out
}

#' @section FacileDGEResult:
#' Given a FacileDGEResult, we can re-materialize the Bioconductor assay
#' container used within the differential testing pipeline used from [fdge()].
#'
#' @export
biocbox.FacileDGEResult <- function(x, ...) {
  res <- biocbox(x[["model_def"]], x[["assay_name"]], x[["method"]],
                 filter = dge(x)$feature_id, ...)
  res
}

#' @noRd
#' @export
dge <- function(x, ...) {
  assert_class(x, "FacileDGEResult")
  multiGSEA::logFC(x[["result"]])
}

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

#' @noRd
#' @export
vizualize.FacileTtestDGEResult <- function(x, type = c("volcano", "feature"),
                                           max_padj = 0.01, min_abs_logFC = 1,
                                           feature = NULL, event_source = "A",
                                           ...) {
  type <- match.arg(type)

}

#' @section Interacting with results:
#'
#' The `report` function will create an htmlwidget which can be explored by
#' the analyst or dropped into an Rmarkdown report.
#'
#' `report(result, "dge", max_padj = 0.05, min_abs_logFC = 1)` will create a
#' side-by-side volcano and datatable for differential expression results.
#'
#' @export
#' @rdname fdge
report.FacileTtestDGEResult <- function(x, result = "dge",
                                        max_padj = 0.01, min_abs_logFC = 1,
                                        event_source = "A", ...) {
  # result <- match.arg(result, c("dge", multiGSEA::resultNames(result$gsea)))
  result <- match.arg(result, "dge")

}

# Helpers ======================================================================

.dge_methods <- c("voom", "qlf", "trended", "limma")

#' A table of assay_type,dge_method combination parameters
#'
#' This table is used to match the assay_type,dge_method combination with
#' the appropriate bioc container class `bioc_class`, and default edgeR/limma
#' model fitting params.
#'
#' @noRd
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
