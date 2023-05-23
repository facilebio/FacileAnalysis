# Things fail when there is multi-modal data stored in the faciledataset, and
# the assay you run the analysis on isn't available on all of the samples you
# sent into the analysis.
#
# This is because `initialized.FacilePcaAnalysisResult` tests to see that the
# number of rows in the output from `tidy(fpca(samples)) == nrow(samples)`, and
# if some rows in `samples` don't have data from the assay under test, then we
# get hosed
library(FacileData)

# efds <- FacileData::exampleFacileDataSet()e
efds <- FacileDataSet("~/workspace/facilebio/data/FacileKPMPDataSet/")
user <- Sys.getenv("USER")
options(facile.log.level.fshine = "trace")
devtools::load_all(".")


# debug(FacileData:::biocbox.facile_frame)
# debug(FacileData:::.fetch_assay_data)

# All in one module ============================================================
shiny::shinyApp(
  ui = shiny::fluidPage(
     FacileShine::filteredReactiveFacileDataStoreUI("rfds"),
    # shiny::tags$h2("fdgeAnalysis"),
    # fdgeAnalysisUI("analysis"),
    shiny::tags$h2("fpcaAnalysis"),
    fpcaAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(FacileShine::filteredReactiveFacileDataStore, "rfds",
                       path = reactive(efds$parent.dir),
                       user = "lianoglou")
    # analysis <- callModule(fdgeAnalysis, "analysis", rfds)
    analysis <- callModule(fpcaAnalysis, "analysis", rfds)
  }
)

if (FALSE) {
  asamples <- samples(efds) |> collect(n = Inf)
  kpca <- fpca(asamples, "scrnaseq")
}