library(FacileData)
library(FacileShine)

fds <- FacileData::exampleFacileDataSet()
user <- Sys.getenv("USER")

devtools::load_all(".")

shiny::shinyApp(
  ui = shiny::fluidPage(
    shiny::wellPanel(filteredReactiveFacileDataStoreUI("ds")),
    shiny::tags$h2("fdgeModelDef"),
    fdgeModelDefUI("model"),
    NULL),
  server = function(input, output) {
    rfds <- callModule(filteredReactiveFacileDataStore, "ds", fds, user = user)
    model <- callModule(fdgeModelDef, "model", rfds)
  }
)
