#' Builds (simple) design and contrast matrices for use with [fdge()]
#'
#' This simplifies the design and contrast building process by allowing for
#' simple model definitions that are, essentially, functions of a single
#' covariate. More elaborate models can be analysed, but the user is left to
#' define the design, coef / contrast to test manually and pass those into
#' [fdge()].
#'
#' Note: actually a (likely) small modification of this can have it support the
#' "ratio of ratios" model setup.
#'
#' @section Missing Covariates:
#' Given the "ragged" nature of sample annotations in a FacileDataStore, some
#' samples may have NA's as their values for the covariates under test. In this
#' case. In this case, if `on_missing` is set to "error", an error will be
#' thrown, otherwise a message will be set in the `warning` list element.
#'
#' The samples that the differential expression should be run on will be
#' enumerated by the `(dataset,sample_id)` pair in the `result$covariates`
#' tibble.
#' 
#' @section Alignment with assay data:
#' This builds a linear model by working with the covariates that are defined
#' over the samples. This does not ask which assay will be used downstream in
#' combination with this linear model to run the fit and test. It is the
#' responsibility of the downstream users/functions of this linear model to
#' ensure that the samples defined in the linear model have data from the
#' assay that the actual measurements/data is coming from.
#'
#' @export
#'
#' @param x a dataset
#' @param covariate the name of the "main effect" sample_covariate we are
#'   performing a contrast against.
#' @param numer character vector defining the covariate/groups that
#'   make up the numerator
#' @param denom character vector defining the covariate/groups that
#'   make up the denominator
#' @param batch character vector defining the covariate/groups to
#'   use as batch effects
#' @param block a string that names the covariate to use for the blocking
#'   factor in a random effects model.
#' @param on_missing when a covariate level is missing (NA) for a sample, the
#'   setting of this parameter (default `"warn"`) will dictate the behavior
#'   of this funciton. When `"warning"`, a warning will be raised, and the
#'   message will be stored in the `$warning` element of the resul. Otherwise,
#'   when `"error"`. See the "Missing Covariates" section for more information.
#'
#' @return a list with:
#'   * `$test`: "ttest" or "anova"
#'   * `$covariates`: the pData over the samples (datset,sample_id, ...)
#'   * `$design`: the design matrix (always 0-intercept)
#'   * `$contrast`: the contrast vector that defines the comparison asked for
#'   * `$messages`: A character vector of messages generated
#'   * `$warnings`: A character vector of warnings generated
#'   * `$errors`: A character vector of errors generated
#'
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#'
#' # Look for tumor vs normal differences, controling for stage and sex
#' model_info <- efds |>
#'   FacileData::filter_samples(indication == "BLCA") |>
#'   flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
#'           batch = "sex")
#' m2 <- efds |>
#'   FacileData::filter_samples(indication == "BLCA") |>
#'   flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
#'           batch = c("sex", "stage"))
#'
#' # stageIV vs stageII & stageIII
#' m3 <- efds |>
#'   FacileData::filter_samples(indication == "BLCA", sample_type == "tumor") |>
#'   flm_def(covariate = "stage", numer = "IV", denom = c("II", "III"),
#'           batch = "sex")
#'
#' # Incomplete ttest to help with custom contrast vector
#' mi <- efds |>
#'   FacileData::filter_samples(indication == "BLCA", sample_type == "tumor") |>
#'   flm_def(covariate = "stage", batch = "sex", contrast. = "help")
#'
#' # ANOVA across stage in BLCA, control for sex
#' m3 <- efds |>
#'   FacileData::filter_samples(indication == "BLCA") |>
#'   flm_def(covariate = "stage", batch = "sex")
flm_def <- function(x, covariate, numer = NULL, denom = NULL,
                    batch = NULL, block = NULL,
                    on_missing = c("warning", "error"), ...) {
  UseMethod("flm_def", x)
}

#' @export
#' @rdname flm_def
#' @method flm_def data.frame
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts nonEstimable
#' @importFrom FacileViz unselected
#' @importFrom stringr str_detect
#'
#' @section data.frame:
#' The `*.data.frame` function definition assumes that `x` is a data.frame of
#' samples (dataset,sample_id) and the covariates defined on these samples
#' (ie. all the other columns of `x`) contain a superset of the variable names
#' used in the construction of the design matrix for the model definition.
#'
#' @param contrast. A custom contrast vector can be passed in for extra tricky
#'   comparisons that we haven't figured out how to put a GUI in front of.
flm_def.data.frame <- function(x, covariate, numer = NULL, denom = NULL,
                               batch = NULL, block = NULL,
                               on_missing = c("warning", "error"), ...,
                               contrast. = NULL,
                               .fds = NULL) {
  on_missing <- match.arg(on_missing)
  assert_subset(c("dataset", "sample_id"), colnames(x))
  assert_choice(covariate, setdiff(colnames(x), c("sample_id")))
  assert_categorical(x[[covariate]])
  if (!is.null(contrast.)) assert_string(contrast.)

  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(batch)) batch <- NULL
  if (unselected(block)) block <- NULL

  assert_subset(batch, setdiff(colnames(x), c("sample_id")))
  if (!is.null(.fds)) {
    assert_facile_data_store(.fds)
  }

  all_test_levels <- local({
    vals <- x[[covariate]]
    if (is.factor(vals)) levels(droplevels(vals)) else sort(unique(vals))
  })
  if (length(all_test_levels) == 2L && is.null(numer) && is.null(denom)) {
    # If this is specified as an ANOVA (no numer or denom), we still run it as a
    # t-test
    numer <- all_test_levels[2L]
    denom <- all_test_levels[1L]
  }

  test_levels <- assert_subset(c(numer, denom), all_test_levels)

  messages <- character()
  warnings <- character()
  errors <- character()

  out <- list(
    # Standard FacileAnalysisResult things
    fds = .fds)
  class(out) <- c("FacileLinearModelDefinition", "FacileAnalysisResult")
  clazz <- "IncompleteModelDefintion"
  on.exit({
    out[["messages"]] <- messages
    out[["warnings"]] <- warnings
    out[["errors"]] <- errors
    class(out) <- c(clazz, class(out))
    return(out)
  })

  if (is.null(test_levels) && is.null(contrast.)) {
    test_type <- "anova"
  } else {
    test_type <- "ttest"
  }

  x <- distinct(x, dataset, sample_id, .keep_all = TRUE)
  # Build the design matrix ----------------------------------------------------
  req.cols <- c("dataset",  "sample_id", covariate, batch, block)
  incomplete <- !complete.cases(select(x, !!req.cols))
  if (any(incomplete)) {
    msg <- paste(sum(incomplete), "samples with NA's in required covariates")
    if (on_missing == "error") stop(msg)
    warning(msg)
    warnings <- c(warnings, paste("Removed", msg))
    x <- x[!incomplete,]
  }

  # Avoids adding all-0 columns to the design matrix, which come from levels
  # of the tested covariate that do not appear in our sample space
  x <- droplevels(x)

  dformula <- paste0(
    "~ ",
    if (test_type == "anova") NULL else "0 + ",
    covariate)

  if (!is.null(batch)) {
    dformula <- paste(dformula, "+", paste(batch, collapse = " + "))
  }

  design <- model.matrix(formula(dformula), data = x)
  # safe_covname <- make.names(covariate)
  # test_covs <- grep(paste0("^", make.names(safe_covname)), colnames(design))

  # The indices of the columns of the design matrix that come from the covariate
  # under test
  test_covs <- grep(paste0("^", covariate), colnames(design))

  # remove the covariate prefix from the colnames of the matrix
  colnames(design) <- sub(paste0("^", covariate), "", colnames(design))

  # Once we strip the covariate prefix from the main effect, we may again
  # introduce invalid colnames, ie. if the covariate was genotype, and one
  # value is 5xFAD, then if we strip off covariate, the column will just
  # Protect against non-valid column names. Do not mangle the first (Intercept)
  # which will be there if this is an ANOVA
  if (test_type == "anova") {
    colnames(design)[-1L] <- make.names(colnames(design)[-1L])
  } else {
    colnames(design) <- make.names(colnames(design))
  }

  rownames(design) <- paste(x$dataset, x$sample_id, sep = "__")

  non_estimable <- nonEstimable(design)
  if (!is.null(non_estimable)) {
    err <- paste("Design matrix is not full rank.",
                 "Cannot estimate these covariates:",
                 paste(non_estimable, collapse = ","),
                 "\n")
    if (is.null(batch)) {
      err <- glue(err, "This is a catastrophic error, please contact ",
                  "lianoglou@dnli.com for help")
    } else {
      check <- sapply(batch, function(fcov) any(grepl(fcov, non_estimable)))
      if (length(check)) {
        err <- glue(err,
                    "Try removing one of these covariates from the model: ",
                    paste(batch[check], collapse = ","))
      } else {
        err <- glue(err, "There must be a problem in sample annotation, ",
                    "please contact lianoglou@dnli.com")
      }
    }
    errors <- c(errors, err)
  }

  if (test_type == "ttest") {
    dup.terms <- intersect(numer, denom)
    if (length(dup.terms)) {
      warnings <- c(
        warnings,
        glue("Warning: the same term(s) appear in both numerator and ",
             "denominator (", paste(dup.terms, collapse = ","), ")")
      )
    }
    if (is.null(contrast.)) {
      if (unselected(numer) || unselected(denom)) {
        errors <- c(
          errors,
          "T-test requires both numerator and denominator to be specified")
      }
    } else {
      # User passed in a custom contrast -- do we need to do anything?
    }
  }

  numer. <- denom. <- contrast_string <- NULL
  if (length(errors)) {
    clazz <- "FacileFailedModelDefinition"
    coef <- NULL
    contrast <- NULL
  } else if (test_type == "anova") {
    clazz <- "FacileAnovaModelDefinition"
    coef <- 2L:max(test_covs)
    contrast <- NULL
    contrast_string <- NULL
    numer. <- denom. <- NULL
  } else {
    coef <- NULL
    if (is.null(contrast.)) {
      numer. <- paste(make.names(numer), collapse = " + ")
      if (length(numer) > 1L) {
        numer. <- sprintf("( %s ) / %d", numer., length(numer))
      }
      denom. <- paste(make.names(denom), collapse = " + ")
      if (length(denom) > 1L) {
        denom. <- sprintf("( %s ) / %d", denom., length(denom))
      }
      contrast_string <- sprintf("%s - %s", numer., denom.)
    } else {
      numer. <- "__inferred__"
      denom. <- "__inferred__"
      contrast_string <- contrast.
    }
    clazz <- "FacileTtestModelDefinition"
    # Are we testing an interaction term, ie. you can test the differences in
    # the treatment effect of a drug1 in hek cells vs its effect in hela cells:
    # (treatment_hek - ctrl_hek) - (treatment_hela - ctrl_hela)
    # Interaction regex, ie. (treat1_ - ctrl) - (treat2 - ctrl)
    iregex <- "\\(.+-.+\\)\\W*-\\W*\\(.+-.+\\)"

    if (str_detect(contrast_string, iregex)) {
      clazz <- c("FacileInteractionTestModelDefinition", clazz)
    }
    contrast <- makeContrasts(contrasts = contrast_string, levels = design)[,1L]
  }

  out[["params"]] <- list(
    covariate = covariate,
    numer = numer,
    denom = denom,
    batch = batch,
    block = block,
    contrast. = contrast.)
  # out[["test_type"]] <- test_type
  out[["covariates"]] <- x
  out[["numer"]] <- numer.
  out[["denom"]] <- denom.
  out[["design_formula"]] <- dformula
  out[["design"]] <- design
  out[["test_covs"]] <- test_covs
  out[["coef"]] <- coef
  out[["contrast"]] <- contrast
  out[["contrast_string"]] <- contrast_string
  out
}

#' @export
#' @rdname flm_def
#' @importFrom FacileViz unselected
flm_def.tbl <- function(x, covariate, numer = NULL, denom = NULL,
                        batch = NULL, block = NULL,
                        on_missing = c("warning", "error"), ...) {
  x <- collect(x, n = Inf)
  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(batch)) batch <- NULL

  flm_def.data.frame(x, covariate = covariate, numer = numer,
                     denom = denom, batch = batch, block = block,
                     on_missing = on_missing, ...)
}

#' @section facile_frame:
#' When we define a model off of a facile_frame, we expect this to look like
#' a wide covariate table. This defines the samples we will build a model on
#' in its (datset, sample_id) columns, as well as any covaraites defined on
#' these samples.
#'
#' If there are covariates used in the `covariate` or `batch` parameters that
#' are not found in `colnames(x)`, we will attempt to retrieve them from the
#' FacileDataStore `fds(x)`. If they cannot be found, this function will raise
#' an error.
#'
#' @export
#' @importFrom FacileViz unselected
#' @rdname flm_def
flm_def.facile_frame <- function(x, covariate, numer = NULL, denom = NULL,
                                 batch = NULL, block = NULL,
                                 on_missing = c("warning", "error"), ...,
                                 custom_key = NULL) {
  .fds <- assert_class(fds(x), "FacileDataStore")
  assert_sample_subset(x)

  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(batch)) batch <- NULL
  if (unselected(block)) block <- NULL

  if (!is.null(block)) assert_string(block)

  # Retrieve any covariates from the FacileDataStore that are not present
  # in the facile_frame `x`
  required.covs <- setdiff(
    c(covariate, batch, block),
    c("dataset", "sample_id"))
  assert_character(required.covs)

  fetch.covs <- setdiff(required.covs, colnames(x))
  if (length(fetch.covs)) {
    x <- with_sample_covariates(x, fetch.covs, custom_key = custom_key,
                                .fds = .fds)
  }

  x <- collect(x, n = Inf)

  out <- flm_def.data.frame(x, covariate = covariate, numer = numer,
                            denom = denom, batch = batch, block = block,
                            on_missing = on_missing, .fds = .fds, ...)
  out
}

#' @export
#' @rdname flm_def
#' @importFrom FacileViz unselected
flm_def.FacileDataStore <- function(x, covariate, numer = NULL, denom = NULL,
                                    batch = NULL, block = NULL,
                                    on_missing = c("warning", "error"),
                                    ..., samples = NULL, custom_key = NULL) {
  if (is.null(samples)) samples <- samples(x)
  samples <- collect(samples, n = Inf)

  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(batch)) batch <- NULL
  if (unselected(block)) block <- NULL

  flm_def(samples, covariate = covariate, numer = numer, denom = denom,
          batch = batch, block = block, on_missing = on_missing,
          custom_key = custom_key, ...)
}

#' @export
#' @rdname flm_def
#' @importFrom FacileShine active_samples
#' @importFrom FacileViz unselected
flm_def.ReactiveFacileDataStore <- function(x, covariate, numer = NULL,
                                            denom = NULL, batch = NULL,
                                            block = NULL,
                                            on_missing = c("warning", "error"),
                                            ..., samples = active_samples(x),
                                            custom_key = user(x)) {
  samples <- collect(samples, n = Inf)
  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(batch)) batch <- NULL
  if (unselected(block)) block <- NULL

  flm_def(samples, covariate = covariate, numer = numer, denom = denom,
          batch = batch, block = block, on_missing = on_missing,
          custom_key = custom_key, ...)
}

# Accessor Functions ===========================================================

#' @noRd
#' @export
samples.FacileTtestModelDefinition <- function(x, tested_only = FALSE, ...) {
  out <- x[["covariates"]]
  if (tested_only) {
    cov <- param(x, "covariate")
    num <- param(x, "numer")
    den <- param(x, "denom")
    out <- filter(out, .data[[cov]] %in% c(num, den))
  }
  out
}


#' @noRd
#' @export
samples.FacileLinearModelDefinition <- function(x, ...) {
  x[["covariates"]]
}

#' @noRd
#' @export
design.FacileLinearModelDefinition <- function(x, ...) {
  x[["design"]]
}

#' @noRd
#' @export
result.FacileLinearModelDefinition <- function(x, ...) {
  x[["design"]]
}

#' The `formula()` function here just returns the formula string for the
#' design matrix
#' @noRd
#' @importFrom stats formula
#' @export
formula.FacileLinearModelDefinition <- function(x, ..., as.string = TRUE) {
  out <- x[["design_formula"]]
  if (!as.string) {
    stop("What should we be doing with this if not return the string?")
  }
  out
}

#' @noRd
.with_warnings <- function(x, fresult, ...) {
  assert_class(x, "FacileAnalysisResultStatus")
  assert_class(fresult, "FacileAnalysisResult")
  w <- fresult[["warnings"]]
  if (length(w)) {
    cc <- class(x)
    x <- paste0(x, "\n", "**Warnings**\n", paste(w, collapse = "\n"))
    class(x) <- cc
  }
  x
}

#' @noRd
#' @export
status.FacileFailedModelDefinition <- function(x, type = "message", ...) {
  out <- paste(x$errors, collapse = "\n")
  class(out) <- c("FacileAnalysisStatusError", "FacileAnalysisResultStatus",
                  "character")
  out
}

#' @noRd
#' @export
status.FacileAnovaModelDefinition <- function(x, type = "message",
                                              with_warnings = TRUE, ...) {
  nsamples <- nrow(samples(x))
  out <- sprintf("ANOVA model defined across '%s' covariate on %d samples",
                 param(x, "covariate"), nsamples)
  class(out) <- c("FacileAnalysisStatusSuccess", "FacileAnalysisResultStatus",
                  "character")
  if (with_warnings) out <- .with_warnings(out, x)
  out
}

#' @noRd
#' @export
status.FacileTtestModelDefinition <- function(x, type = "message",
                                              with_warnings = TRUE, ...) {
  nsamples <- nrow(samples(x))
  out <- sprintf("T-test model [%s] defined on %d samples",
                 x$contrast_string, nsamples)
  class(out) <- c("FacileAnalysisStatusSuccess", "FacileAnalysisResultStatus",
                  "character")
  if (with_warnings) out <- .with_warnings(out, x)
  out
}

#' @noRd
#' @export
status.FacileInteractionTestModelDefinition <- function(x, type = "message",
                                                        with_warnings = TRUE,
                                                        ...) {
  out <- sprintf("Interaction model defined on %d samples", nsamples)
  class(out) <- c("FacileAnalysisResultStatus", "character")
  if (with_warnings) out <- .with_warnings(out, x)
  out
}

# Printing =====================================================================

#' @noRd
#' @export
print.FacileLinearModelDefinition <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
format.FacileLinearModelDefinition <- function(x, ...) {
  if (is_ttest(x)) {
    testing <- sprintf("Testing contrast over `%s` covariate:\n  %s",
                       param(x, "covariate"),
                       x[["contrast_string"]])
  } else {
    testing <- sprintf("Running ANOVA over the levels of the `%s` coefficient",
                       param(x, "covariate"))
  }
  out <- paste(
    "===========================================================\n",
    sprintf("%s\n", class(x)[1L]),
    "-----------------------------------------------------------\n",
    "Design: ", formula(x), "\n",
    testing, "\n",
    "===========================================================\n",
    sep = "")
  out
}

