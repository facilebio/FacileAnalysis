context("PCA (fpca)")

test_that("fpca simply doesn't explode", {
  m <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(m) <- letters[1:10]
  pdat <- as.data.frame(FacileAnalysis:::example_aes_data_table(10, n.cats = 5))
  rownames(pdat) <- colnames(m)

  res <- fpca(m, col_covariates = pdat, use_irlba = FALSE)
  expect_equal(nrow(result(res)), ncol(m))
  expect_setequal(res[["feature_stats"]][["PC"]], colnames(res[["rotation"]]))
})

