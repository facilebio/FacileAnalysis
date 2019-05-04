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
#' @return A `FacileDGEModelDefinition` object, the output from
#'   [fdge_model_def()].
fdgeModelDefRun <- function(input, output, session, rfds, ...,
                            debug = FALSE, .reactive = TRUE) {
  isolate. <- if (.reactive) base::identity else shiny::isolate

  active.samples <- reactive({
    req(initialized(rfds))
    isolate.(active_samples(rfds))
  })

  testcov <- callModule(categoricalSampleCovariateSelect, "testcov",
                        rfds, ..., .with_none = FALSE,
                        .reactive = .reactive)

  numer <- callModule(categoricalSampleCovariateLevels, "numer",
                      rfds, testcov, .reactive = .reactive)

  denom <- callModule(categoricalSampleCovariateLevels, "denom",
                      rfds, testcov, .reactive = .reactive)

  fixedcov <- callModule(categoricalSampleCovariateSelect, "fixedcov",
                         rfds, ..., .with_none = FALSE, .reactive = .reactive)

  model <- reactive({
    # TODO: everytime a change is made, this fires twice
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
    req(!partial)

    out <- fdge_model_def(samples., testcov., numer = numer., denom = denom.,
                          fixed = fixed.)
    out
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
fdgeModelDefRunUI <- function(id, ..., debug = FALSE) {
  ns <- NS(id)

  out <- tagList(
    fluidRow(
      column(
        3,
        categoricalSampleCovariateSelectUI(
          ns("testcov"),
          label = "Grouping Covariate",
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
          label = "Fixed",
          multiple = TRUE)))
  )

  if (debug) {
    out <- tagList(
      out,
      shiny::verbatimTextOutput(ns("debug"), placeholder = TRUE))
  }

  out
}
