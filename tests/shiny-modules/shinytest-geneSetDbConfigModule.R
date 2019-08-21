library(FacileData)
library(shiny)

devtools::load_all(".")
options(facile.log.level.fshine = "trace")

efds <- FacileData::exampleFacileDataSet()
gdb.base <- multiGSEA::getMSigGeneSetDb("h", "human", "entrez")
ttest.res <- FacileData::exampleFacileDataSet() %>%
  FacileData::filter_samples(indication == "CRC") %>%
  fdge_model_def(covariate = "sample_type",
                 numer = "tumor", denom = "normal", batch = "sex") %>%
  fdge(method = "voom")

modres <- FacileAnalysis::frunGadget(
  geneSetDbConfig, geneSetDbConfigUI, x = ttest.res, aresult = ttest.res,
  gdb = gdb.base, retval = "gdb")

