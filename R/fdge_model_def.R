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
#' @export
#'
#' @param x a dataset
#' @param covariate the name of the "main effect" sample_covariate we are
#'   performing a contrast against.
#' @param numer character vector defining the covariate/groups that
#'   make up the numerator
#' @param denom character vector defining the covariate/groups that
#'   make up the denominator
#' @param fixed character vector defining the covariate/groups to
#'   use as fixed effects
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
#' model_info <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor",
#'                  denom = "normal",
#'                  fixed = "sex")
#' m2 <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor",
#'                  denom = "normal",
#'                  fixed = c("sex", "stage"))
#'
#' # stageIV vs stageII & stageIII
#' m3 <- efds %>%
#'   filter_samples(indication == "BLCA", sample_type == "tumor") %>%
#'   fdge_model_def(covariate = "stage", numer = "IV", denom = c("II", "III"),
#'                  fixed = "sex")
#'
#' # Incomplete ttest to help with custom contrast vector
#' mi <- efds %>%
#'   filter_samples(indication == "BLCA", sample_type == "tumor") %>%
#'   fdge_model_def(covariate = "stage", fixed = "sex", contrast. = "help")
#'
#' # ANOVA across stage in BLCA, control for sex
#' m3 <- efds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "stage", fixed = "sex")
fdge_model_def <- function(x, covariate, numer = NULL, denom = NULL,
                           fixed = NULL, on_missing = c("warning", "error"),
                           ...) {
  UseMethod("fdge_model_def", x)
}

#' @export
#' @rdname fdge_model_def
#' @method fdge_model_def data.frame
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts nonEstimable
#' @importFrom FacileShine unselected
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
fdge_model_def.data.frame <- function(x, covariate, numer = NULL, denom = NULL,
                                      fixed = NULL,
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
  if (unselected(fixed)) fixed <- NULL

  assert_subset(fixed, setdiff(colnames(x), c("sample_id")))
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
  class(out) <- c("FacileDGEModelDefinition", "FacileAnalysisResult")
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

  # Build the design matrix ----------------------------------------------------
  req.cols <- c("dataset",  "sample_id", covariate, fixed)
  xx <- x[, req.cols]
  incomplete <- !complete.cases(xx)
  if (any(incomplete)) {
    msg <- paste(sum(incomplete), "samples with NA's in required covariates")
    if (on_missing == "error") stop(msg)
    warning(msg)
    warnings <- c(warnings, paste("Removed", msg))
    xx <- xx[!incomplete,]
  }

  dformula <- paste(
    "~",
    if (test_type == "anova") "" else "0 +",
    covariate)

  if (!is.null(fixed)) {
    dformula <- paste(dformula, "+", paste(fixed, collapse = " + "))
  }

  design <- model.matrix(formula(dformula), data = xx)
  # Protect against non-valid column names
  colnames(design) <- make.names(colnames(design))

  safe_covname <- make.names(covariate)
  test_covs <- grep(paste0("^", make.names(safe_covname)), colnames(design))
  colnames(design) <- sub(paste0("^", safe_covname), "", colnames(design))
  rownames(design) <- paste(xx$dataset, xx$sample_id, sep = "__")

  non_estimable <- nonEstimable(design)
  if (!is.null(non_estimable)) {
    err <- paste("Design matrix is not full rank.",
                 "Cannot estimate these covariates:",
                 paste(non_estimable, collapse = ","),
                 "\n")
    if (is.null(fixed)) {
      err <- glue(err, "This is a catastrophic error, please contact ",
                  "lianoglou@dnli.com for help")
    } else {
      check <- sapply(fixed, function(fcov) any(grepl(fcov, non_estimable)))
      if (length(check)) {
        err <- glue(err,
                    "Try removing one of these covariates from the model: ",
                    paste(fixed[check], collapse = ","))
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
        glue("Warning: the same term(s) apper in both numerator and ",
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
    clazz <- "FacileTtestDGEModelDefinition"
    # Are we testing an interaction term, ie. you can test the differences in
    # the treatment effect of a drug1 in hek cells vs its effect in hela cells:
    # (treatment_hek - ctrl_hek) - (treatment_hela - ctrl_hela)
    # Interaction regex, ie. (treat1_ - ctrl) - (treat2 - ctrl)
    iregex <- "\\(.+-.+\\)\\W*-\\W*\\(.+-.+\\)"
    if (str_detect(contrast_string, iregex)) {
      clazz <- c("FacileInteractionTestDGEModelDefinition", clazz)
    }
    contrast <- makeContrasts(contrasts = contrast_string, levels = design)[,1L]
  }

  out[["params"]] <- list(
    covariate = covariate,
    numer = numer,
    denom = denom,
    fixed = fixed,
    contrast. = contrast.)
  out[["test_type"]] <- test_type
  out[["covariates"]] <- xx
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
#' @rdname fdge_model_def
#' @importFrom FacileShine unselected
fdge_model_def.tbl <- function(x, covariate, numer = NULL, denom = NULL,
                          fixed = NULL, on_missing = c("warning", "error"),
                          ...) {
  x <- collect(x, n = Inf)
  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(fixed)) fixed <- NULL

  fdge_model_def.data.frame(x, covariate = covariate, numer = numer,
                            denom = denom, fixed = fixed,
                            on_missing = on_missing, ...)
}

#' @section facile_frame:
#' When we define a model off of a facile_frame, we expect this to look like
#' a wide covariate table. This defines the samples we will build a model on
#' in its (datset, sample_id) columns, as well as any covaraites defined on
#' these samples.
#'
#' If there are covariates used in the `covariate` or `fixed` parameters that
#' are not found in `colnames(x)`, we will attempt to retrieve them from the
#' FacileDataStore `fds(x)`. If they cannot be found, this function will raise
#' an error.
#'
#' @export
#' @importFrom FacileShine unselected
#' @rdname fdge_model_def
fdge_model_def.facile_frame <- function(x, covariate, numer = NULL,
                                        denom = NULL,
                                        fixed = NULL,
                                        on_missing = c("warning", "error"), ...,
                                        custom_key = NULL) {
  .fds <- assert_class(fds(x), "FacileDataStore")
  assert_sample_subset(x)

  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(fixed)) fixed <- NULL

  # Retrieve any covariates from the FacileDataStore that are not present
  # in the facile_frame `x`
  required.covs <- setdiff(c(covariate, fixed), c("dataset", "sample_id"))
  assert_character(required.covs)

  fetch.covs <- setdiff(required.covs, colnames(x))
  if (length(fetch.covs)) {
    x <- with_sample_covariates(x, fetch.covs, custom_key = custom_key,
                                .fds = .fds)
  }

  x <- collect(x, n = Inf)

  out <- fdge_model_def.data.frame(x, covariate = covariate, numer = numer,
                                   denom = denom, fixed = fixed,
                                   on_missing = on_missing, .fds = .fds, ...)
  out
}

#' @export
#' @rdname fdge_model_def
#' @importFrom FacileShine unselected
fdge_model_def.FacileDataStore <- function(x, covariate, numer = NULL,
                                           denom = NULL, fixed = NULL,
                                           on_missing = c("warning", "error"),
                                           ...,
                                           samples = NULL,
                                           custom_key = NULL) {
  if (is.null(samples)) samples <- samples(x)
  samples <- collect(samples, n = Inf)

  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(fixed)) fixed <- NULL

  fdge_model_def(samples, covariate = covariate, numer = numer, denom = denom,
                 fixed = fixed, on_missing = on_missing,
                 custom_key = custom_key, ...)
}

#' @export
#' @rdname fdge_model_def
#' @importFrom FacileShine active_samples unselected
fdge_model_def.ReactiveFacileDataStore <- function(x, covariate, numer = NULL,
                                                   denom = NULL, fixed = NULL,
                                                   on_missing = c("warning", "error"),
                                                   ...,
                                                   samples = active_samples(x),
                                                   custom_key = user(x)) {
  samples <- collect(samples, n = Inf)
  if (unselected(numer)) numer <- NULL
  if (unselected(denom)) denom <- NULL
  if (unselected(fixed)) fixed <- NULL

  fdge_model_def(samples, covariate = covariate, numer = numer, denom = denom,
                 fixed = fixed, on_missing = on_missing,
                 custom_key = custom_key, ...)
}

# Accessor Functions ===========================================================

#' @noRd
#' @export
samples.FacileDGEModelDefinition <- function(x, ...) {
  x[["covariates"]]
}

#' @noRd
#' @export
design.FacileDGEModelDefinition <- function(x, ...) {
  x[["design"]]
}

# Printing =====================================================================\

#' @noRd
#' @export
print.FacileDGEModelDefinition <- function(x, ...) {
  cat(format(x, ...), "\n")
}

#' @noRd
#' @export
format.FacileDGEModelDefinition <- function(x, ...) {
  testtype <- x$test_type
  if (testtype == "ttest") {
    testing <- "contrast"
    thetest <- paste(names(x$contrast), x$contrast,
                     sep = ": ", collapse = " || ")
  } else {
    testing <- "coefficient"
    thetest <- paste(colnames(x$design)[x$coef], collapse = "|")
  }
  out <- paste(
    "===========================================================\n",
    sprintf("%s\n", class(x)[1L]),
    "-----------------------------------------------------------\n",
    "Design: ", x[["design_formula"]], "\n",
    sprintf("Testing `%s` (%s):\n    %s\n", x[["covariate"]],
            testing, thetest),
    "===========================================================\n",
    sep = "")
  out
}

