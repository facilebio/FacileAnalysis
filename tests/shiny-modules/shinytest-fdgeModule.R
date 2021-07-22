library(FacileShine)

efds <- FacileData::exampleFacileDataSet()
user <- Sys.getenv("USER")
options(facile.log.level.fshine = "trace")
devtools::load_all(".")

# All in one module ============================================================
shiny::shinyApp(
  ui = shiny::fluidPage(
    # reactiveFacileDataStoreUI("rfds"),
    filteredReactiveFacileDataStoreUI("rfds"),
    shiny::tags$h2("fdgeAnalysis"),
    fdgeAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(filteredReactiveFacileDataStore, "rfds",
                       path = reactive(efds$parent.dir),
                       user = "lianoglou")
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)
  }
)
