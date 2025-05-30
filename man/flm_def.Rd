% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flm_def-base.R
\name{flm_def}
\alias{flm_def}
\alias{flm_def.data.frame}
\alias{flm_def.tbl}
\alias{flm_def.facile_frame}
\alias{flm_def.FacileDataStore}
\title{Builds (simple) design and contrast matrices for use with \code{\link[=fdge]{fdge()}}}
\usage{
flm_def(
  x,
  covariate,
  numer = NULL,
  denom = NULL,
  batch = NULL,
  block = NULL,
  on_missing = c("warning", "error"),
  ...,
  metadata = list()
)

\method{flm_def}{data.frame}(
  x,
  covariate,
  numer = NULL,
  denom = NULL,
  batch = NULL,
  block = NULL,
  on_missing = c("warning", "error"),
  ...,
  metadata = list(),
  contrast. = NULL,
  .fds = NULL
)

\method{flm_def}{tbl}(
  x,
  covariate,
  numer = NULL,
  denom = NULL,
  batch = NULL,
  block = NULL,
  on_missing = c("warning", "error"),
  ...,
  metadata = list()
)

\method{flm_def}{facile_frame}(
  x,
  covariate,
  numer = NULL,
  denom = NULL,
  batch = NULL,
  block = NULL,
  on_missing = c("warning", "error"),
  ...,
  metadata = list(),
  custom_key = NULL
)

\method{flm_def}{FacileDataStore}(
  x,
  covariate,
  numer = NULL,
  denom = NULL,
  batch = NULL,
  block = NULL,
  on_missing = c("warning", "error"),
  ...,
  metadata = list(),
  samples = NULL,
  custom_key = NULL
)
}
\arguments{
\item{x}{a dataset}

\item{covariate}{the name of the "main effect" sample_covariate we are
performing a contrast against.}

\item{numer}{character vector defining the covariate/groups that
make up the numerator}

\item{denom}{character vector defining the covariate/groups that
make up the denominator}

\item{batch}{character vector defining the covariate/groups to
use as batch effects}

\item{block}{a string that names the covariate to use for the blocking
factor in a random effects model.}

\item{on_missing}{when a covariate level is missing (NA) for a sample, the
setting of this parameter (default \code{"warn"}) will dictate the behavior
of this funciton. When \code{"warning"}, a warning will be raised, and the
message will be stored in the \verb{$warning} element of the resul. Otherwise,
when \code{"error"}. See the "Missing Covariates" section for more information.}

\item{contrast.}{A custom contrast vector can be passed in for extra tricky
comparisons that we haven't figured out how to put a GUI in front of.}
}
\value{
a list with:
\itemize{
\item \verb{$test}: "ttest" or "anova"
\item \verb{$covariates}: the pData over the samples (datset,sample_id, ...)
\item \verb{$design}: the design matrix (always 0-intercept)
\item \verb{$contrast}: the contrast vector that defines the comparison asked for
\item \verb{$messages}: A character vector of messages generated
\item \verb{$warnings}: A character vector of warnings generated
\item \verb{$errors}: A character vector of errors generated
}
}
\description{
This simplifies the design and contrast building process by allowing for
simple model definitions that are, essentially, functions of a single
covariate. More elaborate models can be analysed, but the user is left to
define the design, coef / contrast to test manually and pass those into
\code{\link[=fdge]{fdge()}}.
}
\details{
Note: actually a (likely) small modification of this can have it support the
"ratio of ratios" model setup.
}
\section{Missing Covariates}{

Given the "ragged" nature of sample annotations in a FacileDataStore, some
samples may have NA's as their values for the covariates under test. In this
case. In this case, if \code{on_missing} is set to "error", an error will be
thrown, otherwise a message will be set in the \code{warning} list element.

The samples that the differential expression should be run on will be
enumerated by the \verb{(dataset,sample_id)} pair in the \code{result$covariates}
tibble.
}

\section{Alignment with assay data}{

This builds a linear model by working with the covariates that are defined
over the samples. This does not ask which assay will be used downstream in
combination with this linear model to run the fit and test. It is the
responsibility of the downstream users/functions of this linear model to
ensure that the samples defined in the linear model have data from the
assay that the actual measurements/data is coming from.
}

\section{data.frame}{

The \verb{*.data.frame} function definition assumes that \code{x} is a data.frame of
samples (dataset,sample_id) and the covariates defined on these samples
(ie. all the other columns of \code{x}) contain a superset of the variable names
used in the construction of the design matrix for the model definition.
}

\section{facile_frame}{

When we define a model off of a facile_frame, we expect this to look like
a wide covariate table. This defines the samples we will build a model on
in its (datset, sample_id) columns, as well as any covaraites defined on
these samples.

If there are covariates used in the \code{covariate} or \code{batch} parameters that
are not found in \code{colnames(x)}, we will attempt to retrieve them from the
FacileDataStore \code{fds(x)}. If they cannot be found, this function will raise
an error.
}

\examples{
efds <- FacileData::exampleFacileDataSet()

# Look for tumor vs normal differences, controling for stage and sex
model_info <- efds |>
  FacileData::filter_samples(indication == "BLCA") |>
  flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
          batch = "sex")
m2 <- efds |>
  FacileData::filter_samples(indication == "BLCA") |>
  flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
          batch = c("sex", "stage"))

# stageIV vs stageII & stageIII
m3 <- efds |>
  FacileData::filter_samples(indication == "BLCA", sample_type == "tumor") |>
  flm_def(covariate = "stage", numer = "IV", denom = c("II", "III"),
          batch = "sex")

# Incomplete ttest to help with custom contrast vector
mi <- efds |>
  FacileData::filter_samples(indication == "BLCA", sample_type == "tumor") |>
  flm_def(covariate = "stage", batch = "sex", contrast. = "help")

# ANOVA across stage in BLCA, control for sex
m3 <- efds |>
  FacileData::filter_samples(indication == "BLCA") |>
  flm_def(covariate = "stage", batch = "sex")
}
