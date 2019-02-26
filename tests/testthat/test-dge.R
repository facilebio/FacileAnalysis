context("Differential Gene Expression and GSEA")

FDS <- FacileData::exampleFacileDataSet()

test_that("fdge_model_def supports simple t-test specification", {
  fsamples <- filter_samples(FDS, indication == "BLCA")

  mdef <- fdge_model_def(fsamples, covariate = "sample_type",
                         numer = "tumor", denom = "normal",
                         fixed = "sex")

})

test_that("fdge_model_def removes samples with NA in covariates", {
  # Test for differences in "subtype_molecular" samples in all BLCA vs BLCA
  # tumor samples. These two should produce the same stuff, but the former
  # should have a warning about removing samples since the normal samples
  # aren't classified with a subtype.

  # No warnings: BLCA,tumor samples should all have stage information.
  stumor <- filter_samples(FDS, indication == "BLCA", sample_type == "tumor")
  mod.tumor <- fdge_model_def(stumor, covariate = "subtype_molecular",
                              numer = "luminal", denom = "basal",
                              fixed = "sex")

  # This will emit a warning since numor samples don't have stage covariates
  sall <- filter_samples(FDS, indication == "BLCA")
  sall <- with_sample_covariates(sall, "subtype_molecular")
  warn.regex <- paste(sum(is.na(sall$subtype_molecular)), "samples with NA")

  mod.all <- expect_warning({
    fdge_model_def(sall, covariate = "subtype_molecular",
                   numer = "luminal", denom = "basal",
                   fixed = "sex")
  }, warn.regex)
  in.warnings <- sapply(mod.all$warnings, function(w) grepl(warn.regex, w))
  expect_logical(in.warnings, min.len = 1L,
                 info = "warning stored in result$warnings")
  expect_true(any(in.warnings))

  # FacileDGEModel spec is the same after NA samples are removed
  expect_equal(mod.all$covariates, mod.tumor$covariates)
  expect_equal(mod.all$design, mod.tumor$design)
})

test_that("model_def supports simple ANOVA specification", {

})
