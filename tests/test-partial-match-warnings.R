# There is a warning of a partial match of `coef` to `coefficients` somewhere
# when the result of a ttest are presented in the app

library(testthat)
library(dplyr)
library(FacileDenaliDataSets)
# library(FacileAnalysis)
devtools::load_all(".")
mfds <- FacileDenaliDataSet("mouse")

opts <- options(warnPartialMatchArgs = TRUE,
                warnPartialMatchDollar = TRUE,
                warnPartialMatchAttr = TRUE)

xs <- mfds %>%
  filter_samples(study_name == "DST-150",
                 library_strategy == "clonetech-lowinput-3prime") %>%
  with_sample_covariates()

dge.eq <- xs %>%
  flm_def("group", "APP_SAA_Hom__Microglia", "WT__Microglia", batch = "sex") %>%
  fdge(method = "edgeR-qlf")

viz(dge.eq, "ENSMUSG00000027523")

dge.vm <- dge.eq %>%
  model() %>%
  fdge(method = "voom")

viz(dge.vm, "ENSMUSG00000027523")
