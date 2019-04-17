library(FacileData)
library(FacileShine)

efds <- FacileData::exampleFacileDataSet()
dfds <- FacileData::FacileDataSet("/Users/lianoglou/workspace/data/FacileData/denali/FacileDenaliMouseDataSet/")
fds <- dfds

user <- Sys.getenv("USER")

devtools::load_all(".")

# All in one module ============================================================
full <- shiny::shinyApp(
  ui = shiny::fluidPage(
    shiny::wellPanel(filteredReactiveFacileDataStoreUI("ds")),
    shiny::tags$h2("fdgeAnalysis"),
    fdgeAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(filteredReactiveFacileDataStore, "ds", fds, user = user)
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)
  }
)

full
