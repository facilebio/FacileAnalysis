context("Model definition for Gene Expression and GSEA (fdge_model_def)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("fdge_model_def supports simple t-test specification", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor", denom = "normal",
                   fixed = "sex")
  expect_is(mdef, "FacileTtestDGEModelDefinition")
  expect_equal(mdef$contrast, c(normal = -1, tumor = 1, sexf = 0))

  # Flip numerator and denominator
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "normal", denom = "tumor",
                   fixed = "sex")
  expect_is(mdef, "FacileTtestDGEModelDefinition")
  expect_equal(mdef$contrast, c(normal = 1, tumor = -1, sexf = 0))
})

test_that("Partial t-test spec is not allowed (no numer or denom)", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "normal", denom = NULL,
                   fixed = "sex")
  expect_is(mdef, "FacileFailedModelDefinition")
  expect_true(length(mdef$errors) == 1L)
})

test_that("fdge_model_def supports ANOVA specification", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "stage", fixed = "sex")
  expect_is(mdef, "FacileAnovaModelDefinition")
  expect_equal(colnames(mdef$design)[1], "(Intercept)")
  expect_equal(mdef$coef, match(c("II", "III", "IV"), colnames(mdef$design)))
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
  n.na <- sum(is.na(sall$subtype_molecular))
  warn.regex <- paste(n.na, "samples with NA")

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
  expect_is(mod.all, "FacileTtestDGEModelDefinition")
  expect_is(mod.tumor, "FacileTtestDGEModelDefinition")
})

test_that("fdge_model_def supports retrieving test covaraites on the fly", {
  # sample descriptor with all required covariates for model building
  ff0 <- FDS %>%
    filter_samples(indication == "BLCA", sample_type == "tumor") %>%
    with_sample_covariates(c("stage", "sex"))
  mref <- fdge_model_def(ff0, covariate = "stage", fixed = "sex")

  # sample descriptor without covariates under test
  ff1 <- mutate(ff0, sex = NULL)
  mtest <- fdge_model_def(ff1, covariate = "stage", fixed = "sex")

  expect_is(mref, "FacileAnovaModelDefinition")
  expect_is(mtest, "FacileAnovaModelDefinition")

  expect_matrix(mref$design)
  expect_equal(mtest$design, mref$design)

  expect_numeric(mref$coef)
  expect_equal(mtest$coef, mref$coef)
})

test_that("Errors gracefully with duplicate entries in numer and denom", {
  samples <- filter_samples(FDS, indication == "CRC")
  good.model <- expect_warning({
    samples %>%
      fdge_model_def(covariate = "subtype_crc_cms",
                     numer = c("CMS1", "CMS2"),
                     denom = c("CMS3", "CMS4"))
  }, "NA")

  bad1 <- expect_warning({
    samples %>%
      fdge_model_def(covariate = "subtype_crc_cms",
                     numer = "CMS1",
                     denom = c("CMS1", "CMS2", "CMS3"))
  }, "NA.*required covariates")

  bad1 <- expect_warning({
    samples %>%
      fdge_model_def(covariate = "subtype_crc_cms",
                     numer = "CMS1", denom = "CMS1")
  }, "NA.*required covariates")
})


test_that("fdge_model_def errors on non-fullrank matrices", {
  samples <- filter_samples(FDS, indication == "CRC")

  good.model <- expect_warning({
    samples %>%
      fdge_model_def(covariate = "subtype_crc_cms",
                     numer = c("CMS1", "CMS2"),
                     denom = c("CMS3", "CMS4"))
  }, "NA")
  expect_class(good.model, "FacileTtestDGEModelDefinition")

  # adding `fixed = "sex"` makes this not full rank
  bad.model <- expect_warning({
    samples %>%
      fdge_model_def(covariate = "subtype_crc_cms",
                     numer = c("CMS1", "CMS2"),
                     denom = c("CMS3", "CMS4"),
                     fixed = "sex")
  }, "NA")
  expect_class(bad.model, "FacileFailedModelDefinition")
  expect_string(bad.model$errors, pattern = "full rank")
  expect_string(bad.model$errors, pattern = "removing.*model:.*sex$")
})

test_that("invalid R variable named covariate levels are safe", {
  samples <- filter_samples(FDS, indication == "CRC")
  model <- expect_warning({
    fdge_model_def(samples, covariate = "subtype_microsatellite_instability")
  }, "NA.*required covariates")
})
