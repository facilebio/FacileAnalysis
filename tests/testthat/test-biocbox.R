context("Building bioconductor assay containers with biocbox (biocbox)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Correct bioc container built from model def and method spec", {
  mdef <- FDS |>
    filter_samples(indication == "BLCA") |>
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal")

  # DGEList for quasi-likelihood teseting
  y <- biocbox(mdef, "rnaseq", method = "edgeR-qlf")
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
  expect_set_equal(rownames(vm), rownames(y))

  # EList without weights for limma-trend analysis
  e <- biocbox(mdef, "rnaseq", method = "limma-trend")
  expect_class(e, "EList")
  expect_null(e$weights)
  expect_matrix(e$E, nrows = nrow(e), ncols = ncol(e))
  expect_set_equal(rownames(e), rownames(y))

  # the expression matrices returned from voom and limma-trended are the same
  cors <- cor(e$E, vm$E, method = "spearman")
  expect_true(all(diag(cors) > 0.999))
})

test_that("Various filter* combinations to biocbox work", {
  mdef <- FDS |>
    filter_samples(indication == "BLCA") |>
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal")

  all.features <- features(FDS, assay_name = "rnaseq")

  # ELists generated for limma-trend (skip the voom step)

  # No filtering takes place when filter = NULL
  full.box <- biocbox(mdef, "rnaseq", method = "limma-trend", filter = FALSE)
  expect_equal(nrow(full.box), nrow(all.features))
  expect_set_equal(rownames(full.box), all.features$feature_id)

  # Default low-expression filtering happens when filter isn't specified
  filtered.box <- biocbox(mdef, "rnaseq", method = "limma-trend")
  expect_true(nrow(filtered.box) < nrow(full.box))
  expect_subset(rownames(filtered.box), rownames(full.box))

  # Filtering works on a restricted universe.
  # Here we restrict the universe to genes with 5 letter names.
  # A more common usecase might be to filter the features based on
  # meta == "protein_coding"
  fives.all <- filter(all.features, nchar(name) == 5)

  # Low expression filtering happens within restricted universe
  fives.box <- biocbox(mdef, "rnaseq", method = "limma-trend",
                       filter_universe = fives.all)
  expect_true(all(nchar(fives.box$genes$symbol) == 5))
  expect_subset(fives.box$genes$feature_id, fives.all$feature_id)

  # Specifying filter_require rescues lowly-expressed genes
  add.req <- setdiff(fives.all$feature_id, rownames(fives.box))
  add.req <- sample(add.req, 10)
  expect_false(any(add.req %in% rownames(fives.box)))


  rescued.box <- biocbox(mdef, "rnaseq", method = "limma-trend",
                         filter_universe = fives.all, filter_require = add.req)
  expect_true(nrow(rescued.box) == nrow(fives.box) + length(add.req))
  expect_set_equal(rownames(rescued.box), c(rownames(fives.box), add.req))
})

test_that("biocbox(..., features = something) only ever returns `something`", {
  mdef <- FDS |>
    filter_samples(indication == "BLCA") |>
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal")

  all.features <- features(FDS, assay_name = "rnaseq")
  some.features <- sample_n(all.features, 10)

  box1 <- biocbox(mdef, "rnaseq", method = "limma-trend",
                  filter = TRUE, features = some.features)
  expect_set_equal(rownames(box1), some.features$feature_id)

  box2 <- biocbox(mdef, "rnaseq", method = "limma-trend",
                  filter = "default", features = some.features)
  expect_set_equal(rownames(box2), some.features$feature_id)
})
