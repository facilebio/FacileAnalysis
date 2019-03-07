library(FacileData)
library(FacileShine)

fds <- FacileData::exampleFacileDataSet()
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
    shiny::tags$h2("fdgeModelDef"),
    fdgeModelDefUI("model"),
    shiny::tags$h2("fdgeRun"),
    fdgeRunUI("dge"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(filteredReactiveFacileDataStore, "ds", fds, user = user)
    model <- callModule(fdgeModelDef, "model", rfds)
    dge <- callModule(fdgeRun, "dge", rfds, model)
  }
)
