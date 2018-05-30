expect_color_map <- function(map, values, info = "unknown color map test") {
  uvals <- unique(values)
  expect_equal(length(map), length(uvals), info = info)
  expect_true(setequal(names(map), uvals), info = info)
}

example_aes_data_table <- function(n = 20, n.cats = 3, seed = 123) {
  if (is.numeric(seed)) set.seed(seed[1L])
  x <- tibble(
    category = sample(head(letters, n.cats), n, replace = n.cats < n),
    subcat = sample(c("sub_1", "sub_2"), n, replace = n > 2),
    score = rnorm(n))
  x
}
