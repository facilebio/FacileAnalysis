context("Building bioconductor assay containers with biocbox (biocbox)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Correct bioc container built from model def and method spec", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  # DGEList for quasi-likelihood teseting
  bbox <- biocbox(mdef, "rnaseq", method = "edgeR-qlf")
  y <- result(bbox)
  expect_class(y, "DGEList")

  expect_subset("sample_type", colnames(y$samples))

  # expect dispersions estimated
  expect_number(y$common.dispersion)
  expect_numeric(y$trended.dispersion, len = nrow(y))
  expect_numeric(y$tagwise.dispersion, len = nrow(y))


  # EList with $weights matrix for voom
  bbox <- biocbox(mdef, "rnaseq", method = "voom")
  vm <- result(bbox)
  expect_class(vm, "EList")
  expect_matrix(vm$weights, nrows = nrow(vm), ncols = ncol(vm))
  expect_equal(y$design, vm$design)

  # EList without weights for limma-trend analysis
  bbox <- biocbox(mdef, "rnaseq", method = "limma-trend")
  e <- result(bbox)
  expect_class(e, "EList")
  expect_null(e$weights)
  expect_matrix(e$E, nrows = nrow(e), ncols = ncol(e))

  # the expression matrices returned from voom and limma-trended are the same
  cors <- cor(e$E, vm$E, method = "spearman")
  expect_true(all(diag(cors) > 0.999))
})
