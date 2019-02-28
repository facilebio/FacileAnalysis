context("Facile PCA (fpca)")

test_that("fpca works (this isn't really teseting anything)", {
  m <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(m) <- letters[1:10]
  pdat <- as.data.frame(FacileAnalysis:::example_aes_data_table(10, n.cats = 5))
  rownames(pdat) <- colnames(m)

  res <- fpca(m, col_covariates = pdat)
  expect_equal(nrow(res$tidy), ncol(m))
  expect_equal(nrow(res$factor_contrib), nrow(m))
  expect_equal(
    colnames(res$factor_contrib)[1:5],
    c("feature_id", paste0("PC", 1:4)))
})

