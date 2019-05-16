library(FacileDenaliDataSets)
library(FacileShine)

efds <- FacileData::exampleFacileDataSet()
dfds <- FacileDenaliDataSet("mouse", "internal")
fds <- dfds

user <- Sys.getenv("USER")

options(facile.log.level.fshine = "trace")

devtools::load_all(".")

# All in one module ============================================================
shiny::shinyApp(
  ui = shiny::fluidPage(
    reactiveFacileDataStoreUI("rfds"),
    shiny::tags$h2("fdgeAnalysis"),
    fdgeAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    rfds <- ReactiveFacileDataStore(fds, "ds")
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)
  }
)


