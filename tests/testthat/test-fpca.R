# sparse and full samples
if (!exists("ssamples")) ssamples <- FacileData::some_samples(sparse = TRUE)
if (!exists("fsamples")) fsamples <- FacileData::some_samples(sparse = FALSE)

test_that("fpca output has expected things", {
  res <- fpca(fsamples, "scrnaseq", ntop = 100)
  stats <- tidy(res)
  
  expect_setequal(stats$sample_id, fsamples$sample_id)
  
  # Loadings on the PCS are full
  sigs <- tidy(signature(res))
  expect_setequal(sigs$dimension, colnames(res[["rotation"]]))
})

test_that("fpca reports and recovers from running of sparse-sample assays", {
  res <- expect_warning(fpca(ssamples, "scrnaseq", ntop = 100), "no scrnaseq")
  ws <- warnings(res)
  expect_true(any(grepl("no scrnaseq", ws)))
  
  analyzed <- samples(res)
  expect_lt(nrow(analyzed), nrow(ssamples))
  
  scsamples <- ssamples |> 
    has_assay("scrnaseq") |> 
    filter(has_scrnaseq)
  expect_setequal(analyzed$sample_id, scsamples$sample_id)
  
  tres <- tidy(res)  
  expect_setequal(tres$sample_id, analyzed$sample_id)
})

test_that("fpca can run on prespecified features and overrides default filter", {
  use.features <- FacileData::fds(fsamples) |> 
    FacileData::features(assay_name = "scrnaseq") |>
    dplyr::sample_n(20)

  res <- fpca(fsamples, features = use.features)
  expect_setequal(features(res)$feature_id, use.features$feature_id)
  
  res2 <- fpca(fsamples, features = use.features,
               filter = "variance", ntop = 500)
  expect_setequal(features(res2)$feature_id, use.features$feature_id)
})

