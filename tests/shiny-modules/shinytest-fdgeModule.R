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

# Pieces =======================================================================
def_and_res <- shiny::shinyApp(
  ui = shiny::fluidPage(
    shiny::wellPanel(filteredReactiveFacileDataStoreUI("ds")),

    shiny::tags$h2("FacileModule: fdgeModelDef"),
    fdgeModelDefUI("model"),

    shiny::tags$h2("FacileModule: fdgeRun"),
    fdgeRunUI("dge"),

    shiny::tags$h2("FacileModule: fdgeViewResult"),
    fdgeViewResultUI("view"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(filteredReactiveFacileDataStore, "ds", fds, user = user)
    model <- callModule(fdgeModelDef, "model", rfds)
    dge <- callModule(fdgeRun, "dge", rfds, model)
    view <- callModule(fdgeViewResult, "view", rfds, dge)
  }
)
