#' @noRd
#' @export
assay_name.FacileAnalysisResult <- function(x, ...) {
  param(x, "assay_name")
}

#' @noRd
#' @export
assay_type.FacileAnalysisResult <- function(x, ...) {
  assay_info(fds(x), assay_name(x))$assay_type
}

#' @noRd
#' @export
feature_type.FacileAnalysisResult <- function(x, ...) {
  assay_info(fds(x), assay_name(x))$feature_type
}


# Universal FacileAnalysisResult methods =======================================

#' Extract the internal (F)acile(A)nalysis(R)esult (O)bject from a container.
#'
#' `FacileAnalysisResult` object can be stored within a
#' `ReactiveFacileAnalysis`, or a `FacileMultiAnalysisResult` (which is a
#' collection of FacileAnalysisResults returned from a "multi-analysis" view,
#' like [fDgeSeaAnalysis()])
#'
#' Due to the potential hybrid-nature of analysis workflow, we are not always
#' certain that a FacileAnalysisResult used by a module that did not create it
#' (ie. the `dgeres` object in the [fdgeView()] is a
#' `ReactiveFacileAnalysisResult`, or simply an innert `FacileAnalysisResult`.
#'
#' So that we don't have to always check for these things, we call this method
#' so that we get access to the `FacileAnalysisResult` object itself.
#'
#' As an aside, I know this is a bit of a cop-out to doing full-blown OO. I
#' feel like if I were to implement each of the specific `*FacileAnalysisResult`
#' methods for its sister `ReactiveFacile*AnalysisResult`, then this probably
#' wouldn't be necessary ... I think, anyway.
#'
#' @export
#' @return an innert FacileAnalysisResult
faro <- function(x, ...) {
  UseMethod("faro", x)
}

#' @noRd
#' @export
faro.default <- function(x, ...) {
  NULL
}

#' @noRd
#' @export
faro.FacileAnalysisResult <- function(x, ...) {
  x
}

#' Redoes the analysis based on something new (samples, params, etc.)
#' 
#' @export
#' @param x The original FacileAnalysisResult
#' @param ... probably additional arguments to retweak the original analysis
#' @param samples redo the analysis on a subset of these samples
redo <- function(x, ..., samples = NULL) {
  UseMethod("redo", x)
}

#' Reports the status of the AnalysisResult
#'
#' @export
#' @param x A `FacileAnalysisResult`
#' @return a FacileAnalysisResultStatus (message)
status <- function(x, type = "message", ...) {
  UseMethod("status", x)
}

#' consumers of `status` can test if the result of this is a
#' `FacileAnalysisStatus` or not ...
#'
#' @noRd
#' @export
status.default <- function(x, type = "message", ...) {
  "undefined"
}

#' Are two facile objects comparable?
#'
#' Can we compare two analyses? Datastore? etc.
#'
#' @export
comparable <- function(x, y, ...) {
  UseMethod("comparable", x)
}

#' Compares two (or more(?)) FacileAnalysisResults against each other.
#'
#' This was initially motivated by the desire to compare the results of two
#' differential expressio analyses against each other. Extending the idea
#' to comparing GSEA results (using Thomas' enrichmentmap idea) is a natural
#' extension. There may be othres.
#'
#' Not every FacileAnalysisResult can be compared, I guess, in which case it
#' we will throw an error
#'
#' @export
compare <- function(x, y, ...) {
  UseMethod("compare", x)
}

#' Default/core function to compare FacileAnalysisResults by combining the
#' feature-level statistics returned by the ranks induced over each analysis
#' (as long as they are run on the same feature space.)
#'
#' @noRd
#' @export
compare.FacileAnalysisResult <- function(x, y, ...) {
  assert_class(x, "FacileAnalysisResult") # this should be a tautology
  assert_class(y, "FacileAnalysisResult") # this should be a tautology

  # TODO: need to check that x and y are run on same feature_type and induces
  # ranks on them
  # assert_true(feature_type(x) == feature_type(y))
  rx <- ranks(x, ...)
  ry <- ranks(y, ...)
}

#' @export
#' @noRd
metadata <- function(x, ...) {
  UseMethod("metadata", x)
}

#' @export
#' @noRd
metadata.default <- function(x, ...) {
  list()
}

#' @export
#' @noRd
metadata.FacileAnalysisResult <- function(x, ...) {
  x[["metadata"]]
}

#' @export
#' @noRd
result <- function(x, name = "result", ...) {
  UseMethod("result", x)
}

#' @export
#' @noRd
result.FacileAnalysisResult <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' @export
#' @noRd
tidy.FacileAnalysisResult <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}


#' @export
#' @noRd
result.FacileFeatureRanks <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' @export
#' @noRd
tidy.FacileFeatureRanks <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' @export
#' @noRd
result.FacileFeatureSignature <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' @export
#' @noRd
tidy.FacileFeatureSignature <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' Are these results (ranks, signatures) signed?
#'
#' @export
#' @param x a rank or signature result
#' @return logical
signed <- function(x, ...) {
  UseMethod("signed", x)
}

#' @noRd
#' @export
signed.FacileFeatureRanks <- function(x, ...) {
  grepl("Signed$", class(x)[1L], ignore.case = FALSE)
}

#' @noRd
#' @export
signed.FacileFeatureSignature <- function(x, ...) {
  grepl("Signed$", class(x)[1L], ignore.case = FALSE)
}

#' Extract the value of a parameter(s) used in a FacileAnalysis result
#'
#' @export
#' @param x A FacileAnalysisResult
#' @param name the name of the parameter to extract
param <- function(x, name = NULL, ...) {
  UseMethod("param", x)
}

#' @noRd
#' @export
param.FacileAnalysisResult <- function(x, name = NULL, ...) {
  params. <- assert_list(x[["params"]], names = "unique")
  if (is.null(name) || length(name) == 0L) {
    out <- params.
  } else {
    assert_character(name)
    name <- assert_subset(unique(name), names(params.))
    if (length(name) > 1L) {
      out <- params.[name]
    } else {
      out <- params.[[name]]
    }
  }
  out
}

#' @noRd
#' @export
param.BiocBox <- function(x, name = NULL, ...) {
  param.FacileAnalysisResult(x, name, ...)
}

#' Extracts ranks and signatures from a FacileAnalysisResult.
#'
#' (Note: there is a lot of philosophizing going on here).
#' It is often the case that an analysis over a set of features (or samples)
#' induces a ranking over the features (or samples), which is determend by the
#' test performed in the analysis. The `ranks()` and `signatures()` functions
#' returns a ranking induced over the features (or samples) from the analysis.
#'
#' When an analysis imposes ranks, this usually only occurs over only one of
#' "features" or "samples" used in the analysis. In the even that both of these
#' can be ranked, then these functions will accept a `type` argument which you
#' can parameterize with by either `"features"` or `"samples"`.
#'
#' Signatures are essentially a summary extracted from the ranks. This is most
#' often the "topn" ranks returned from the analysis (or a dimension thereof).
#'
#' @section Signed and Unsigned Ranks:
#' Let's consider a differential gene expression (DGE) analysis, where we are
#' testing the differential abundance of a gene across two groups of samples.
#' The result of the analysis can induce both a signed and unsigned ranking on
#' the genes under test.
#'
#' ```
#' dge <- FacileData::exampleFacileDataSet() |>
#'   FacileData::filter_samples(indication == "BLCA") |>
#'   flm_def(covariate = "sample_type",numer = "tumor", denom = "normal",
#'           batch = "sex") |>
#'   fdge()
#' ```
#'
#' Ranking the results of the DGE by ascending *p*-value will provide an
#' **unsigned** ranking on the genes: Alternatively, one could get a **signed**
#' ranking from the result of this test simply by arranging each gene by its
#' log-fold-change. Both approaches are achieved by the code below:
#'
#' ```
#' uranks <- ranks(dge, signed = FALSE)
#' sranks <- ranks(dge, signed = TRUE)
#' ```
#'
#' rank each gene by its *p*-value (ascending).
#'
#' @section Ranks (Differential Expression):
#' Any differential gene expression analysis works over a set of samples in
#' order to find which genes are most differentially abundant between the two
#' conditions that are defined by the sample groups. It's clear to see, here,
#' that this analysis induces a ranking over the *features* (genes). This
#' ranking can be both signed or unsigned.
#'
#'
#' @section Signature (Differential Expression):
#'
#' @section Ranks (PCA):
#'
#' @section Signature (PCA):
#'
#' @section Signatures and Ranks for all things:
#'
#' Think about how the next analysis you implement fits this scenario.
#'
#' Perhaps we can consider an e-/p-/etc- *QTL analysis is on whose features are
#' SNPs, which we can rank by ones that have strong association with the
#' quantitative phenotypes under test.
#'
#' @section Why is this a formalism in the FacileVerse:
#' I feel like being able to generate succinct summaries of an analysis will
#' more easily enable the analyst (via GUI or code) to dive back in and ask
#' another question. That question might be as simple as "how does this version
#' of my question compare to a slightly different version?"
#'
#' @export
#' @param x A `FacileAnalysisResult`
#' @rdname ranks_and_signatures
ranks <- function(x, ...) {
  UseMethod("ranks", x)
}


#' Converts an analysis result into a "signature" type of thing.
#'
#' This is mainly for feature rankings, but could be something else? Returns
#' a table of features that can be piped into GeneSetDb() constructor, should
#' insipiration strike.
#'
#' Most of the `signature.*` functions defined over the result of analysis first
#' passes through a call to `ranks()` first, which then calls the `signature()`
#' method on those ranks, ie. the two segments below will produce the same
#' result.
#'
#' ```r
#' sigs1 <- FacileData::exampleFacileDataSet() |>
#'   FacileData::filter_samples(indication == "CRC") |>
#'   fpca() |>
#'   signature()
#' sigs2 <- FacileData::exampleFacileDataSet() |>
#'   FacileData::filter_samples(indication == "CRC") |>
#'   fpca() |>
#'   ranks() |>
#'   signature()
#' ```
#'
#' @export
#' @param x A `FacileAnalysisResult`.
#' @examples
#' pca.sigs <- FacileData::exampleFacileDataSet() |>
#'   FacileData::filter_samples(indication == "CRC") |>
#'   fpca() |>
#'   signature()
signature <- function(x, ...) {
  UseMethod("signature", x)
}

# Facile API Implementations ===================================================

#' FacileAnalysisResult objects should be able to fetch the FacileDataStore
#' they were created from.
#'
#' @export
#' @noRd
fds.FacileAnalysisResult <- function(x) {
  return(x[["fds"]])
}

#' @noRd
#' @export
organism.FacileAnalysisResult <- function(x, ...) {
  organism(fds(x), ...)
}

# Vizualization and Rmarkdon reporting =========================================

#' Methods to interactively explore and report FacileAnalysisResults.
#'
#' The `vizualize`, `shine`, and `report` triumverate provide the analyst with
#' the tools required to interact and explore the results of a FacileAnalysis.
#'
#' @section Vizualize:
#' The `vizualize` functions generate an analysis-specific interactive
#' htmlwidget for the analysist to explore.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#'
#' @param x A `FacileAnalysisResult` object
#' @param ... passed down to the `x`-specific `vizualize.*`, `report.*`, and
#'   `shine.*` functions.
viz <- function(x, ...) {
  UseMethod("viz", x)
}

#' @section Report:
#' The `report` function produces an object that can be embedded into an
#' Rmarkdown document. The implementation of these functions will likely
#' result in a parameterized call to the respective `vizualize` function.
#'
#' Two `report` functions should be created per FacileAnalysis. One that accepts
#' the `FacileAnalysisResult` object itself, and another that accepts the
#' analysis-specific `FacileAnalysisShine` object, which should parameterize
#' the final `visualize` call so that it can be embedded seamlessly into
#' an Rmarkdown report for "offline" viewing.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#' @aliases report
report <- function(x, ...) {
  UseMethod("report", x)
}

# Statistical Modeling Stuff ===================================================

#' Extract the design matrix from objects used in a linear modeling step
#'
#' @export
design <- function(x, ...) {
  UseMethod("design", x)
}

#' @noRd
#' @export
design.NULL <- function(x, ...) NULL

#' Extract the model used for an analysis
#'
#' For a DGE analysis, this is the flm_def object
#'
#' @export
model <- function(x, ...) {
  UseMethod("model", x)
}

#' @noRd
#' @export
model.NULL <- function(x, ...) NULL

#' Extract the contrast defined.
#' 
#' @export
#' @param x a facileanalysis result / FacileLinearModelDefinition
contrast <- function(x, ...) {
  UseMethod("contrast", x)
}

#' @noRd
#' @export
contrast.default <- function(x, ...) {
  stop("contrast() not defined on object class: ", class(x)[1L])
}

#' @noRd
#' @export
contrast.NULL <- function(x, ...) NULL


# Vizualization and Rmarkdon reporting =========================================

#' Methods to interactively explore and report FacileAnalysisResults.
#'
#' The `vizualize`, `shine`, and `report` triumverate provide the analyst with
#' the tools required to interact and explore the results of a FacileAnalysis.
#'
#' @section Vizualize:
#' The `vizualize` functions generate an analysis-specific interactive
#' htmlwidget for the analysist to explore.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#'
#' @param x A `FacileAnalysisResult` object
#' @param ... passed down to the `x`-specific `vizualize.*`, `report.*`, and
#'   `shine.*` functions.
viz <- function(x, ...) {
  UseMethod("viz", x)
}
