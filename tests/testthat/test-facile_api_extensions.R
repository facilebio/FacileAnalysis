context("Facile API ++")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("with_sample_covs gives precedence to local covariates", {
  xs <- filter_samples(FDS, indication == "BLCA") %>% collect()
  all.covs <- summary(fetch_sample_covariates(xs))

  xs$stage <- sample(letters, nrow(xs), replace = TRUE)

  test1 <- with_sample_covs(xs)
  expect_equal(select(test1, !!colnames(xs)), xs)

  extra <- xs %>%
    mutate(stage = rnorm(nrow(xs)), stuff = sample(letters, nrow(xs)))

  test2 <- with_sample_covs(xs, custom_covariates = extra)
  expect_equal(select(test2, !!colnames(extra)), extra)
})
