context("PCA (fpca)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("fpca simply doesn't explode", {
  res <- FDS |>
    FacileData::filter_samples(indication == "BLCA") |>
    fpca(ntop = 1000, prior.count = 0.5)

  expect_equal(nrow(tidy(res)), nrow(samples(res)))
  expect_setequal(res[["feature_stats"]][["PC"]], colnames(res[["rotation"]]))
})

test_that("fpca can run on prespecified features and overrides default filter", {
  use.features <- FDS |>
    FacileData::features() |>
    dplyr::sample_n(20)

  res <- FDS |>
    FacileData::filter_samples(indication == "BLCA") |>
    fpca(features = use.features, ntop = 1000, prior.count = 0.5)

  features.used <- res |>
    FacileData::features() |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  expect_setequal(features.used$feature_id, use.features$feature_id)
})
