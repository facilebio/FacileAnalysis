context("Building bioconductor assay containers with biocbox (biocbox)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Correct bioc container built from model def and method spec", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  # DGEList for quasi-likelihood teseting
  y <- biocbox(mdef, "rnaseq", method = "qlf")
  expect_class(y, "DGEList")
  expect_subset("sample_type", colnames(y$samples))

  # expect dispersions estimated
  expect_number(y$common.dispersion)
  expect_numeric(y$trended.dispersion, len = nrow(y))
  expect_numeric(y$tagwise.dispersion, len = nrow(y))


  # EList with $weights matrix for voom
  vm <- biocbox(mdef, "rnaseq", method = "voom")
  expect_class(vm, "EList")
  expect_matrix(vm$weights, nrows = nrow(vm), ncols = ncol(vm))
  expect_equal(y$design, vm$design)

  # EList without weights for limma-trend analysis
  e <- biocbox(mdef, "rnaseq", method = "trended")
  expect_class(e, "EList")
  expect_null(e$weights)
  expect_matrix(e$E, nrows = nrow(e), ncols = ncol(e))
  expect_equal(e$E, vm$E)
})
