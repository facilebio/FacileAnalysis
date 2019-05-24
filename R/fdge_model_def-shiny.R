#' A shiny module that generates a linear model definition via `fdge_model_def`.
#'
#' ```
#' model_info <- fds %>%
#'   filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor",
#'                  denom = "normal",
#'                  fixed = "sex")
#' ```
#'
#' @export
#' @importFrom FacileShine
#'   active_samples
#'   categoricalSampleCovariateSelect
#'   categoricalSampleCovariateLevels
#'   initialized
#'   update_exclude
#' @importFrom shiny renderText
#' @importFrom shinyjs toggleElement
#' @importFrom shinyWidgets sendSweetAlert
#' @return A `FacileDGEModelDefinition` object, the output from
#'   [fdge_model_def()].
fdgeModelDefRun <- function(input, output, session, rfds, ...,
                            debug = FALSE, .reactive = TRUE) {
  isolate. <- if (.reactive) base::identity else shiny::isolate

  active.samples <- reactive({
    req(initialized(rfds))
    # isolate.(active_samples(rfds))
    ftrace("Updating active samples")
    active_samples(rfds)
  })

  testcov <- callModule(categoricalSampleCovariateSelect, "testcov",
                        rfds, include1 = FALSE, ..., .with_none = FALSE,
                        .reactive = .reactive)
  # the "fixed" covariate is what I'm calling the extra/batch-level
  # covariates. the entry selected in the testcov is removed from the
  # available elemetns to select from here
  fixedcov <- callModule(categoricalSampleCovariateSelect, "fixedcov",
                         rfds, ..., .with_none = FALSE,
                         .exclude = testcov$covariate,
                         reactive = .reactive)

  numer <- callModule(categoricalSampleCovariateLevels, "numer",
                      rfds, testcov, .reactive = .reactive)

  denom <- callModule(categoricalSampleCovariateLevels, "denom",
                      rfds, testcov, .reactive = .reactive)

  # Make the levels available in the numer and denom covariates
  # mutually exclusive
  # Note: I can't get categoricalSampleCovariateLevels to work
  # smoothly like this when `ignoreNULL = FALSE` is set so that
  # when one selectize drains, its last level is made availalbe
  # to the other select.
  # TODO: Use FacileShine::categoricalSampleCovariateLevels
  #       instead of individual numer and denom selects so that
  #       the empty select releases its last level to the
  #       "select pool"
  # observe({
  #   update_exclude(denom, numer$values)
  #   update_exclude(numer, denom$values)
  # })

  model <- reactive({
    req(initialized(rfds))
    samples. <- active.samples()
    testcov. <- name(testcov)
    req(!unselected(testcov.))

    numer. <- numer$values()
    denom. <- denom$values()
    fixed. <- name(fixedcov)

    # Ensure that either
    #   i. neither numer or denom is filled so that we run an ANOVA;
    #  ii. both are filled for a propper t-test specification
    partial <- xor(unselected(numer.), unselected(denom.))
    all.dups <- !unselected(numer.) && setequal(numer., denom.)

    if (partial || all.dups) {
      out <- NULL
    } else {
      out <- fdge_model_def(samples., testcov., numer = numer., denom = denom.,
                            fixed = fixed.)
    }
    out
  })

  status <- reactive({
    req(initialized(rfds))
    model. <- model()
    if (!is(model., "FacileDGEModelDefinition")) {
      out <- "Uninitialized"
    } else if (is(model., "FacileFailedModelDefinition")) {
      out <- "error"
    } else if (length(model.$warnings) == 0) {
      out <- "warn"
    } else {
      out <- "clean"
    }
    out
  })

  observeEvent(status(), {
    status. <- status()
    # model. <- model()
    # show.mbox <- status. != "error" && is(model, "FacileDGEModelDefinition")
    # toggleElement("messagebox", condition = show.mbox)

    toggleElement("messagebox", condition = status. != "error")
    if (status. == "error") {
      sendSweetAlert(session, "Error building model",
                     text = model()$errors, type = "error")
    }
  })

  output$message <- renderText({
    model. <- model()
    status. <- isolate(status())
    clazz <- class(model.)[1L]
    if (!is.null(model.)) {
      nsamples <- nrow(samples(model.))
    }

    msg <- switch(
      clazz,
      "NULL" = "Undefined model",
      FacileAnovaModelDefinition = {
        sprintf("ANOVA model defined across '%s' covariate on %d samples",
                param(model., "covariate"), nsamples)
      },
      FacileInteractionTestDGEModelDefinition = {
        sprintf("Interaction model defined on %d samples", nsamples)
      },
      FacileTtestDGEModelDefinition = {
        sprintf("T-test model [%s] defined on %d samples",
                model.$contrast_string, nsamples)
      },
      FacileFailedModelDefinition = {
        paste("ERROR:", paste(model.$errors, collapse = "\n"))
      },
      "Unknown model type defined")

    if (status. == "warn") {
      msg <- paste(msg, model.$warnings, sep = "\n")
    }
    msg
  })

  if (debug) {
    output$debug <- shiny::renderText({
      model. <- req(model())
      format(model.)
    })
  }

  vals <- list(
    result = model,
    testcov = testcov,
    numer = numer,
    denom = denom,
    fixedcov = fixedcov,
    .ns = session$ns)

  class(vals) <- c("ShinyDGEModelDefinition", class(vals))
  vals
}

#' @noRd
#' @importFrom FacileShine
#'   categoricalSampleCovariateSelectUI
#'   categoricalSampleCovariateLevelsUI
#' @importFrom shiny textOutput wellPanel
#' @importFrom shinyjs hidden
fdgeModelDefRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  out <- tagList(
    fluidRow(
      column(
        3,
        categoricalSampleCovariateSelectUI(
          ns("testcov"),
          label = "Group to test",
          multiple = FALSE)),
      column(
        3,
        categoricalSampleCovariateLevelsUI(
          ns("numer"),
          label = "Numerator",
          multiple = TRUE)),
      column(
        3,
        categoricalSampleCovariateLevelsUI(
          ns("denom"),
          label = "Denominator",
          multiple = TRUE)),
      column(
        3,
        categoricalSampleCovariateSelectUI(
          ns("fixedcov"),
          label = "Control for",
          multiple = TRUE))),
    hidden(wellPanel(id = ns("messagebox"), textOutput(ns("message"))))
  )

  if (debug) {
    out <- tagList(
      out,
      shiny::verbatimTextOutput(ns("debug"), placeholder = TRUE))
  }

  out
}
