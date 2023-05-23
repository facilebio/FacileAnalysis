efds <- FacileData::exampleFacileDataSet()
user <- Sys.getenv("USER")
options(facile.log.level.fshine = "trace")
devtools::load_all(".")


# debug(FacileData:::biocbox.facile_frame)
# debug(FacileData:::.fetch_assay_data)

# All in one module ============================================================
shiny::shinyApp(
  ui = shiny::fluidPage(
    # reactiveFacileDataStoreUI("rfds"),
    FacileShine::filteredReactiveFacileDataStoreUI("rfds"),
    shiny::tags$h2("fdgeAnalysis"),
    fdgeAnalysisUI("analysis"),
    NULL),

  server = function(input, output) {
    rfds <- callModule(FacileShine::filteredReactiveFacileDataStore, "rfds",
                       path = reactive(efds$parent.dir),
                       user = "lianoglou")
    analysis <- callModule(fdgeAnalysis, "analysis", rfds)
  }
)
