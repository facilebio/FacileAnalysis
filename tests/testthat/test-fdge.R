context("Differential Gene Expression and GSEA")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Simple treatment vs control analysis runs", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  dge_test <- fdge(mdef, method = "qlf", gsea = NULL)

})

test_that("Correct bioc container built from model def and method spec", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  # DGEList for quasi-likelihood teseting
  y <- bioc_container(mdef, "rnaseq", method = "qlf")
  expect_class(y, "DGEList")
  # expect dispersions estimated
  expect_number(y$common.dispersion)
  expect_numeric(y$trended.dispersion, len = nrow(y))
  expect_numeric(y$tagwise.dispersion, len = nrow(y))

  # EList with $weights matrix for voom
  vm <- bioc_container(mdef, "rnaseq", method = "voom")
  expect_class(vm, "EList")
  expect_matrix(vm$weights, nrows = nrow(vm), ncols = ncol(vm))
  expect_equal(y$design, vm$design)

  # EList without weights for limma-trend analysis
  e <- bioc_container(mdef, "rnaseq", method = "trended")
  expect_class(e, "EList")
  expect_null(e$weights)
  expect_matrix(e$E, nrows = nrow(e), ncols = ncol(e))
  expect_equal(e$E, vm$E)
})
