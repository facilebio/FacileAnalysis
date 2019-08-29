#' @noRd
#' @export
#' @importFrom multiGSEA iplot geneSets
viz.FacileFseaAnalysisResult <- function(x, name = NULL, collection = NULL,
                                         ...) {
  mgres <- assert_class(result(x), "MultiGSEAResult")

  # Passing in a geneset name will use multiGSEA::iplot
  if (is.character(name)) {
    stopifnot(length(name) == 1L)
    gs <- geneSets(mgres)
    found <- gs[["name"]] == name
    colxn <- gs[["collection"]][found]
    if (length(colxn) == 0L) {
      stop("Can not find geneset with name: ", name)
    }
    if (length(colxn) > 1L) {
      stop("Multiple genesets found with name '", name, "'. ",
           "Need to specify specify collection: ",
           paste(colxn, collapse = ","))
    }
    out <- list(
      plot = iplot(mgres, colxn, name, ...),
      input_data = NULL,
      params = list(name = name, collection = colxn))
    class(out) <- c("FacileViz")
  }
  out
}

#' @noRd
#' @export
#' @examples
#' gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
#' dge.ttest <- FacileData::exampleFacileDataSet() %>%
#'   FacileData::filter_samples(indication == "CRC") %>%
#'   flm_def(covariate = "sample_type", numer = "tumor", denom = "normal",
#'           batch = "sex") %>%
#'   fdge(method = "voom")
#' gsea.ttest <- ffsea(dge.ttest, gdb, methods = c("cameraPR", "fgsea"))
#'
#' viz(gsea.ttest, name = "HALLMARK_ANGIOGENESIS")
#' if (interactive()) {
#'   shine(gsea.ttest)
#' }
shine.FacileFseaAnalysisResult <- function(x, user = Sys.getenv("USER"),
                                           title = "FSEA Results",
                                           ...) {
  frunGadget(ffseaView, ffseaViewUI, x, aresult = x, title = title, ...)
}
