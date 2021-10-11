#' @noRd
#' @export
viz.FacileFseaAnalysisResult <- function(x, type = c("density", "gsea"),
                                         name = NULL, collection = NULL,
                                         rank_by = NULL, ...) {
  mgres <- assert_class(result(x), "SparrowResult")
  type <- match.arg(type)

  assert_string(name)
  assert_string(collection, null.ok = TRUE)

  if (is.null(rank_by)) {
    rank_by <- param(x, "rank_by")
  }
  assert_string(rank_by)
  if (rank_by != param(x, "rank_by")) {
    warning("Plot generated using different rank statistic than tested")
  }

  # TODO: This needs refactoring to make better use of sparrow::geneSet
  #       retrieval mojo already there .....................................
  # gs <- sparrow::geneSets(mgres)
  # found <- gs[["name"]] == name
  # if (!any(found)) {
  #   stop("No geneset found with name: ", name)
  # }
  # if (!is.null(collection)) {
  #   found <- found & gs[["collection"]] == collection
  # } else {
  #   collection <- gs[["collection"]][found]
  #   if (length(collection) == 0L) {
  #     stop("Can not find geneset with name: ", name)
  #   }
  #   if (length(collection) > 1L) {
  #     stop("Multiple genesets found with name '", name, "'. ",
  #          "Need to specify specify collection: ",
  #          paste(collection, collapse = ","))
  #   }
  # }
  # ........................................................................
  plt <- sparrow::iplot(mgres, name = name, collection = collection,
                        type = type, value = rank_by, ...)

  out <- list(
    plot = plt,
    input_data = NULL,
    params = list(name = name, collection = collection))
  class(out) <- c("FacileViz")
  out
}

#' @noRd
#' @export
#' @examples
#' gdb <- sparrow::getMSigGeneSetDb("h", "human", id.type = "entrez")
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
                                           viewer = "browser", ...) {
  frunGadget(ffseaView, ffseaViewUI, x, aresult = x, title = title,
             viewer = viewer, ...)
}
