#' Builds (simple) design and contrast matrices for use with [fdge()]
#'
#' This simplifies the design and contrast building process by allowing for
#' simple model definitions that are, essentially, functions of a single
#' covariate. More elaborate models can be analysed, but the user is left to
#' define the design, coef / contrast to test manually and pass those into
#' [dge()].
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
#' @param covariate the name of the sample_covariate we are performing a
#'   contrast against
#' @param on_missing when a covariate level is missing (NA) for a sample, the
#'   setting of this parameter (default `"warn"`) will dictate the behavior
#'   of this funciton. When `"warning"`, a warning will be raised, and the
#'   message will be stored in the `$warning` element of the resul. Otherwise,
#'   when `"error"`. See the "Missing Covariates" section for more information.
#' @return a list with:
#'   * `$test`: "ttest" or "anova"
#'   * `$covariates`: the pData over the samples (datset,sample_id, ...)
#'   * `$design`: the design matrix (always 0-intercept)
#'   * `$contrast`: the contrast vector that defines the comparison asked for
#'   * `$messages`: A character vector of messages generated
#'   * `$warnings`: A character vector of warnings generated
#' @examples
#'
#' fds <- FacileData::exampleFacileDataSet()
#'
#' # Look for tumor vs normal differences, controling for stage and sex
#' model_info <- fds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "sample_type", numer = "tumor", denom = "normal",
#'             fixed = c("sex"))
fdge_model_def <- function(x, covariate, numer = NULL, denom = NULL,
                           fixed = NULL, on_missing = c("warning", "error"),
                           ...) {
  UseMethod("fdge_model_def", x)
}

#' @export
#' @rdname fdge_model_def
#' @method fdge_model_def data.frame
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts
#' @section data.frame:
#' The data.frame function definition assumes that `x` is a data.frame of
#' samples (dataset,sample_id) and the covariates defined on these samples
#' that contains a superset of the covariates used in the design and contrast
#' construction.
fdge_model_def.data.frame <- function(x, covariate, numer = NULL, denom = NULL,
                                      fixed = NULL,
                                      on_missing = c("warning", "error"), ...) {
  on_missing <- match.arg(on_missing)
  assert_subset(c("dataset", "sample_id"), colnames(x))
  assert_choice(covariate, setdiff(colnames(x), c("sample_id")))
  assert_categorical(x[[covariate]])
  assert_subset(fixed, setdiff(colnames(x), c("sample_id")))

  all_test_levels <- as.character(unique(x[[covariate]]))
  test_levels <- assert_subset(c(numer, denom), all_test_levels)

  messages <- character()
  warnings <- character()
  errors <- character()

  test <- if (is.null(test_levels)) "anova" else "ttest"

  req.cols <- c("dataset",  "sample_id", covariate, fixed)
  xx <- x[, req.cols]
  incomplete <- !complete.cases(xx)
  if (any(is.na(incomplete))) {
    msg <- paste(sum(incomplete), " samples due to NA's in required covariates")
    if (on_missing == "error") stop(msg)
    warnings <- c(warnings, paste("Removed", msg))
    xx <- xx[!incomplete,]
  }

  dformula <- paste("~ 0 +", covariate)
  if (!is.null(fixed)) {
    dformula <- paste(dformula, "+", paste(fixed, collapse = " + "))
  }

  design <- model.matrix(formula(dformula), data = xx)
  test_covs <- grep(paste0("^", covariate), colnames(design))
  colnames(design) <- sub(paste0("^", covariate), "", colnames(design))

  if (test == "anova") {
    coef <- 2:max(test_covs)
    contrast <- NULL
  } else {
    coef <- NULL
    numer. <- paste(numer, collapse = " + ")
    if (length(numer) > 1L) {
      numer. <- sprintf("(%s) / %d", numer., length(numer))
    }
    denom. <- paste(denom, collapse = " + ")
    if (length(denom.) > 1L) {
      denom. <- sprintf("(%s) / %d", denom., length(denom))
    }
    contrast <- makeContrasts(
      contrasts = sprintf("%s - %s", numer., denom.),
      levels = design)
    contrast <- contrast[, 1L]
  }

  out <- list(
    covariates = xx,
    covariate = covariate,
    fixed = fixed,
    numer = numer,
    denom = denom,
    design_formula = dformula,
    design = design,
    test_covs = test_covs,
    coef = coef,
    contrast = contrast,
    messages = messages,
    warnings = warnings,
    errors = errors)

  class(out) <- "FacileDGEModelDefinition"
  out
}


#' @export
#' @rdname fdge_model_def
fdge_model_def.tbl <- function(x, covariate, numer = NULL, denom = NULL,
                          fixed = NULL, on_missing = c("warning", "error"),
                          ...) {
  x <- collect(x, n = Inf)
  fdge_model_def.data.frame(x, covariate = covariate, numer = numer,
                            denom = denom, fixed = fixed,
                            on_missing = on_missing, ...)
}

#' @export
#' @rdname fdge_model_def
fdge_model_def.FacileDataStore <- function(x, covariate, numer = NULL,
                                           denom = NULL, fixed = NULL,
                                           on_missing = c("warning", "error"),
                                           ...,
                                           samples = active_samples(x),
                                           custom_key = NULL) {
  samples <- collect(samples, n = Inf)
  fdge_model_def(samples, covariate = covariate, numer = numer, denom = denom,
                 fixed = fixed, on_missing = on_missing,
                 custom_key = custom_key, ...)
}

#' @section facile_frame
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
#' @rdname fdge_model_def
fdge_model_def.facile_frame <- function(x, covariate, numer = NULL,
                                        denom = NULL,
                                        fixed = NULL,
                                        on_missing = c("warning", "error"), ...,
                                        custom_key = NULL) {
  .fds <- fds(x) %>% assert_class("FacileDataStore")
  assert_sample_subset(x)

  # Retrieve any covariates from the FacileDataStore that are not present
  # in the facile_frame `x`
  required.covs <- c(covariate, fixed)
  assert_character(required.covs)

  fetch.covs <- setdiff(required.covs, colnames(x))

  if (length(fetch.covs)) {
    x <- with_sample_covariates(x, fetch.covs, custom_key = custom_key,
                                .fds = .fds)
  }

  x <- collect(x, n = Inf)

  out <- fdge_model_def.data.frame(x, covariate = covariate, numer = numer,
                                   denom = denom, fixed = fixed,
                                   on_missing = on_missing, ...)
  out
}
