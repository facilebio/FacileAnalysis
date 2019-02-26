context("Differential Gene Expression and GSEA")

FDS <- FacileData::exampleFacileDataSet()

test_that("fdge_model_def supports simple t-test specification", {
  fsamples <- filter_samples(FDS, indication == "BLCA")

  mdef <- fdge_model_def(fsamples, covariate = "sample_type",
                         numer = "tumor", denom = "normal",
                         fixed = "sex")

})

test_that("model_def supports simple ANOVA specification", {

})
