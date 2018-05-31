context("Facile PCA (fpca)")

test_that("fpca dimension reduction works", {
  m <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(m) <- letters[1:10]
  pdat <- as.data.frame(FacileAnalysis:::example_aes_data_table(10, n.cats = 5))
  rownames(pdat) <- colnames(m)

  res <- fpca(m, col_covariates = pdat)
})

