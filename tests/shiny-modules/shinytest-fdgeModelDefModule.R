library(FacileData)
library(FacileShine)

fds <- FacileData::exampleFacileDataSet()
user <- Sys.getenv("USER")

path <- fds[["parent.dir"]]

devtools::load_all(".")

shiny::shinyApp(
  ui = shiny::fluidPage(
    # shiny::wellPanel(filteredReactiveFacileDataStoreUI("ds")),
    reactiveFacileDataStoreUI("rfds"),
    shiny::tags$h2("fdgeModelDef"),
    fdgeModelDefRunUI("model", debug = TRUE),
    NULL),
  server = function(input, output) {
    # rfds <- callModule(filteredReactiveFacileDataStore, "ds", fds, user = user)
    rfds <- ReactiveFacileDataStore(path, "ds")# , samples = s)
    model <- callModule(fdgeModelDefRun, "model", rfds, debug = TRUE)
  }
)
