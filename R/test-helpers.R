# These functions sometimes might go into tests/testhat/helper-*.R files, but
# I like having them available (and not exported) within the package itself
# in case I want to use them elsewhere.

#' Test for a valid discrete aes map
#'
#' @noRd
#' @param expected name of the expected RColorBrewer palette used in map. If
#'   this is provided, then the colors in map are checked to be the same top n
#'   colors provided by the given map name
expect_daes_map <- function(map, values, expected = NULL,
                            info = "unknown color map test") {
  assert_categorical(values)
  uvals <- if (is.factor(values)) levels(values) else sort(unique(values))
  n.cats <- length(unique(uvals))

  expect_equal(length(map), n.cats, info = info)
  expect_true(setequal(names(map), uvals), info = info)

  if (!is.null(expected)) {
    stopifnot(is.character(expected) || is.integerish(expected))
    if (is.character(expected)) {
      if (is.brewer.map.name(expected)) {
        # percolates a warning if n.cats > maximum color of brewer palette
        expected <- RColorBrewer::brewer.pal(n.cats, expected)
        names(expected) <- head(uvals, length(expected))
      }
    }
    expect_equal(length(expected), n.cats)
    # if (is(map, "AsIs")) {
    #   nms <- names(map)
    #   map <- if (is.categorical(map)) as.character(map) else as.numeric(map)
    #   names(map) <- nms
    # }
    expect_equal(map, expected, info = info)
  }
}

example_aes_data_table <- function(n = 20, n.cats = 3, seed = 123) {
  if (is.numeric(seed)) set.seed(seed[1L])
  x <- tibble(
    category = sample(head(letters, n.cats), n, replace = n.cats < n),
    subcat = sample(c("sub_1", "sub_2"), n, replace = n > 2),
    score = rnorm(n))
  x
}

