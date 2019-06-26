# Universal FacileAnalysisResult methods =======================================

#' Extract the internal (F)acile(A)nalysis(R)esult (O)bject from a
#' ReactiveFacileAnalysis
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

#' @export
#' @return an innert FacileAnalysisResult
faro <- function(x, ...) {
  UseMethod("faro", x)
}

#' @noRd
#' @export
faro.default <- function(x, ...) {
  stop("faro.default is undefined")
}

#' @noRd
#' @export
faro.FacileAnalysisResult <- function(x, ...) {
  x
}

#' @noRd
#' @export
faro.ReactiveFacileAnalysisResult <- function(x, ...) {
  assert_class(x$faro, "reactive")
  x$faro()
}

#' ReactiveFacileAnalysisResultContainer are results from a shiny module that
#' encapsulate a single complete analysis step, such as the fdgeAnalysis, for
#' instance.
#'
#' These `*Analysis` will store the main FacileAnalysisResult in their
#' `"main"` list element.
#'
#' @noRd
#' @export
faro.ReactiveFacileAnalysisResultContainer <- function(x, main = "main", ...) {
  # assert_choice(main, names(x))
  req(test_choice(main, names(x)))
  main. <- x[[main]]
  # assert_class(main., "ReactiveFacileAnalysisResult")
  req(test_class(main., "ReactiveFacileAnalysisResult"))
  faro(main.)
}

#' A shiny module that produces a FacileAnalysisResult as its main
#' result will store it as a reactive() in its result$faro element.
#'
#' Triggering its reactivity will result in the return of that object.
#' If the analysis-result has not be generated yet within that shiny module,
#' then triggering it MUST raise an error, due to the module ensuring
#' req() is in all the right places
#'
#' For instance, the fdgeRun module is responsible for everything being
#' set up properly, and when the Run button is hit, it includes a req()
#' in response to make sure the result is built correctly
#'
#' @noRd
#' @export
#' @importFrom FacileShine initialized
initialized.ReactiveFacileAnalysisResult <- function(x, ...) {
  obj <- try(faro(x), silent = TRUE)
  is(obj, "FacileAnalysisResult") && initialized(obj)
}

#' @noRd
#' @export
initialized.FacileAnalysisResult <- function(x, ...) {
  length(x$errors) == 0L
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

compare.FacileAnalysisResult <- function(x, y, ...) {
  msg <- glue("No `compare` method defined for `{class(x)[1L]}` result type.")
}

# Retrieves the main (or alternative) result from a FacileAnalysis
#
# This should always return a tidy data frame of one sort or another.
#
# FacileAnalysisResults should should all aslo have a top-level result
# accessible via obj[["result"]]
#
#' @importFrom multiGSEA result
#' @export result
NULL

# This is defined in multiGSEA@develop
# result <- function(x, name = "result", ...) {
#   UseMethod("result", x)
# }

#' @export
#' @noRd
result.FacileAnalysisResult <- function(x, name = "result", ...) {
  assert_choice(name, names(x))
  x[[name]]
}

#' @export
#' @noRd
result.ReactiveFacileAnalysisResult <- function(x, name = "result", ...) {
  result(faro(x, ...), name = name, ...)
}

#' Extract the value of a parameter used in a FacileAnalysis result
#'
#' @export
#' @param x A FacileAnalysisResult
#' @param name the name of the parameter to extract
param <- function(x, name, ...) {
  UseMethod("param", x)
}

#' @noRd
#' @export
param.FacileAnalysisResult <- function(x, name, ...) {
  params. <- assert_list(x[["params"]], names = "unique")
  assert_choice(name, names(params.))
  params.[[name]]
}

#' @noRd
#' @export
param.ReactiveFacileAnalysisResult <- function(x, name, ...) {
  assert_choice("result", names(x))
  param(x$result(), name, ...)
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
#' dge <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "BLCA") %>%
#'   fdge_model_def(covariate = "sample_type",
#'                  numer = "tumor", denom = "normal", fixed = "sex") %>%
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
#' sigs1 <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   signature()
#' sigs2 <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
#'   ranks() %>%
#'   signature()
#' ```
#'
#' @export
#' @param x A `FacileAnalysisResult`.
#' @examples
#' pca.sigs <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   fpca() %>%
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

#' @noRd
#' @export
viz.FacileGadgetResult <- function(x, ...) {
  viz(result(x), ...)
}

#' @section Shine:
#' The `shine` functions generate shiny gadgets that provide a more interactive
#' view over a `FacileAnalysisResult`. This empowers the analyst to provide more
#' context around the results, likely by leveraging all of the data available
#' within the FacileDataStore.
#'
#' The respective `shine` functions must return a `FacileAnalysisShine` object
#' invisibly to the caller. These should be able to be past into an overladed
#' `report` function, the result of which can be embedded into an Rmarkdown
#' report. In this way the analyst can embed a feature-reduced version of what
#' was observed in the gadget into an Rmarkdown report.
#'
#' I'm not sure how exactly we can do this, but perhaps this will require some
#' code generation that the analyst can copy and paste paste into the Rmarkdown
#' document. This might simply be a parameterized version of the `report`
#' function call.
#'
#' @export
#' @rdname FacileAnalysisResultViz
#' @aliases shine
shine <- function(x, ...) {
  UseMethod("shine", x)
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

#' @noRd
#' @export
report.FacileGadgetResult <- function(x, ...) {
  report(result(x), ...)
}

# Statistical Modeling Stuff ===================================================

#' Extract the design matrix from objects used in a linear modeling step
#'
#' @export
design <- function(x, ...) {
  UseMethod("design", x)
}

#' Extract the model used for an analysis
#'
#' For a DGE analysis, this is the fdge_model_def object
#'
#' @export
model <- function(x, ...) {
  UseMethod("model", x)
}

# Generic API calls on results from running an analysis through a gadget =======

#' #' @noRd
#' report.FacileGadgetResult <- function(x, ...) {
#'   report(result(x), ...)
#' }
#'
#' #' @noRd
#' viz.FacileGadgetResult <- function(x, ...) {
#'   viz(result(x), ...)
#' }
#'
#' #' @noRd
#' shine.FacileGadgetResult <- function(x, ...) {
#'   shine(result(x), ...)
#' }
#'
#' #' @noRd
#' compare.FacileGadgetResult <- function(x, y, ...) {
#'   if (is(y, "FacileGadgetResult")) y <- result(y)
#'   compare(result(x), y)
#' }
