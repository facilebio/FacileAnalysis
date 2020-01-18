context("Building bioconductor assay containers with biocbox (biocbox)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Correct bioc container built from model def and method spec", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def(covariate = "sample_type",
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
  expect_set_equal(rownames(vm), rownames(y))

  # EList without weights for limma-trend analysis
  bbox <- biocbox(mdef, "rnaseq", method = "limma-trend")
  e <- result(bbox)
  expect_class(e, "EList")
  expect_null(e$weights)
  expect_matrix(e$E, nrows = nrow(e), ncols = ncol(e))
  expect_set_equal(rownames(e), rownames(y))

  # the expression matrices returned from voom and limma-trended are the same
  cors <- cor(e$E, vm$E, method = "spearman")
  expect_true(all(diag(cors) > 0.999))
})

test_that("Various filter* combinations to biocbox work", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal")

  all.features <- features(FDS, assay_name = "rnaseq")

  # ELists generated for limma-trend (skip the voom step)

  # No filtering takes place when filter = NULL
  f.all <- mdef %>%
    biocbox("rnaseq", method = "limma-trend", filter = NULL) %>%
    features()
  expect_equal(nrow(f.all), nrow(all.features))
  expect_set_equal(f.all$feature_id, all.features$feature_id)

  # Default low-expression filtering happens when filter isn't specified
  f.default <- mdef %>%
    biocbox("rnaseq", method = "limma-trend") %>%
    features()
  expect_true(nrow(f.default) < nrow(f.all))
  expect_true(all(f.default$feature_id %in% f.all$feature_id))

  # Filtering works on a restricted universe.
  # Here we restrict the universe to genes with 5 letter names.
  # A more common usecase might be to filter the features based on
  # meta == "protein_coding"
  fives.all <- filter(all.features, nchar(name) == 5)

  # Low expression filtering happens within restricted universe
  fives.default <- mdef %>%
    biocbox("rnaseq", method = "limma-trend", filter_universe = fives.all) %>%
    features()
  expect_true(all(nchar(fives.default$symbol) == 5))
  expect_true(all(fives.default$feature_id %in% f.default$feature_id))

  # Specifying filter_require rescues lowly-expressed genes
  add.req <- setdiff(fives.all$feature_id, fives.default$feature_id)
  add.req <- sample(add.req, 10)
  expect_false(any(add.req %in% fives.default$feature_id))

  fives.with.req <- mdef %>%
    biocbox("rnaseq", method = "limma-trend", filter_universe = fives.all,
            filter_require = add.req) %>%
    features()
  expect_set_equal(
    fives.with.req$feature_id,
    c(fives.default$feature_id, add.req),)

  # filter = feature descriptor does that and just that
  take.these <- head(fives.all, 10)
  only.these <- mdef %>%
    biocbox("rnaseq", method = "limma-trend", filter = take.these) %>%
    features()
  expect_set_equal(only.these$feature_id, take.these$feature_id)
})
