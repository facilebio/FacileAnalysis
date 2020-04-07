# library(FacileDenaliDataSets)
# dfds <- FacileDenaliDataSet("mouse", "internal")
# fds <- dfds

library(FacileShine)

efds <- FacileData::exampleFacileDataSet()

user <- Sys.getenv("USER")

options(facile.log.level.fshine = "trace")

devtools::load_all(".")

# All in one module ============================================================
shiny::shinyApp(
  ui = shiny::fluidPage(
    shinyjs::useShinyjs(),
    # reactiveFacileDataStoreUI("rfds"),
    FacileShine:::filteredReactiveFacileDataStoreUI("rfds"),
    shiny::tags$h2("fdgeAnalysis"),
    fdgeAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    # rfds <- ReactiveFacileDataStore(efds, "ds")
    rfds <- callModule(filteredReactiveFacileDataStore, "rfds",
                       reactive(efds$parent.dir))
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)
  }
)


