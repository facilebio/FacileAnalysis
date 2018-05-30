context("Aesthetics")

test_that("create_color_map maps RColorBrewer color palettes to levels", {
  # This tests the simple case where the number of unique elements is fewer
  # than the number of colors in a given palette.
  vals <- head(letters, 5L)
  nvals <- length(vals)

  pal.names <- c("Set1", "Set2", "Set3")
  for (pal in pal.names) {
    expected <- RColorBrewer::brewer.pal(nvals, pal)
    names(expected) <- vals
    map <- create_color_map(vals, pal)
    expect_daes_map(map, vals, expected, info = paste("Small n,", pal))
    # this is a test for the expect_daes_map function
    expect_daes_map(map, vals, pal, info = paste("Small n,", pal))
  }
})

test_that("create_color_map recycle colors when necessary", {
  # When the number of unique levels is more than the colors in the color map
  # we simply expect that:
  # 1. Every element still gets a color
  # 2. Colors will be re-used among different elements
  cols <- FacileAnalysis:::mucho.colors()
  ncols <- length(cols)

  vals <- head(letters, ncols + 5)
  map <- create_color_map(vals, cols)
  dup.cols <- duplicated(map)

  expect_daes_map(map, vals, info = "more elements than colors")

  expect_true( !any(head(dup.cols,  ncols)), info = "First color appearance" )
  expect_true(  all(tail(dup.cols, -ncols)), info = "Duplicated colors" )
})

# aes_map setting and retrieval ================================================

test_that("aes_map() returns empty list on object with no mappings", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)
  y <- rnorm(10)

  amap.x <- aes_map(x)
  expect_is(amap.x, "list")
  expect_length(amap.x, 0)

  amap.y <- aes_map(y)
  expect_is(amap.y, "list")
  expect_length(amap.y, 0)
})

test_that("aes maps can be set individually", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)
  aes_map(x, "color") <- create_color_map(x$category)

  amap <- aes_map(x)
  expect_daes_map(amap$color, x$category)
})

test_that("multiple aes maps can be set by a list", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)
  aes_map(x) <- list(color = create_color_map(x$category), size = "all sizes")

  amap <- aes_map(x)
  expect_daes_map(amap$color, x$category)
  expect_equal(amap$size, "all sizes")
})

test_that("list aes_map setter replaces previously defined aesthetic maps", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)
  aes_map(x, "color") <- create_color_map(x$category, "Set1")
  expect_daes_map(aes_map(x, "color"), x[["category"]], "Set1")
  # expect_equal(
  #   unname(aes_map(x, "color")),
  #   RColorBrewer::brewer.pal(n.cats, "Set1"))


  expect_warning({
    aes_map(x) <- list(color = create_color_map(x$category, "Set3"),
                       size = "all sizes")
  }, "Replacing aesthetic.*color")
  expect_daes_map(aes_map(x, "color"), x[["category"]], "Set3")
  # expect_equal(
  #   unname(aes_map(x, "color")),
  #   RColorBrewer::brewer.pal(n.cats, "Set3"))
  expect_equal(aes_map(x, "size"), "all sizes")
})

# with_* =======================================================================

test_that(".aes_varval_colnames helper function is kosher", {
  n.obs <- 20
  n.cats <- 3
  dat <- example_aes_data_table(n.obs, n.cats)

  test.me <- c(".color_aes.", ".color_aes....", ".color_aes")
  for (vname in test.me) {
    vv <- .aes_varval_colnames(vname)
    expect_equal(vv[["variable"]], ".color_aes.variable", info = vname)
    expect_equal(vv[["value"]], ".color_aes.value", info = vname)
  }
})

test_that("with_color adds color columns to data.frames", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)

  x <- with_color(x, aesthetic = "category", "Set3")
  xc.map <- aes_map(x, "color")
  expect_daes_map(xc.map, x[["category"]], "Set3")

  y <- expect_warning(with_color(x, aesthetic = "category", "Set1"))
  yc.map <- aes_map(y, "color")
  expect_daes_map(yc.map, y[["category"]], "Set1")

  expect_equal(length(xc.map), length(yc.map))
  expect_equal(names(xc.map), names(yc.map))
  expect_false(all(xc.map == yc.map))
})

test_that("with_shape adds color columns to data.frames", {
  n.obs <- 20
  n.cats <- 3
  x <- example_aes_data_table(n.obs, n.cats)

  vv <- .aes_varval_colnames(".shape_aes.")
  x <- with_shape(x, aesthetic = "category")

  # confirm shapes are added
  expect_is(x[[vv$variable]], "character")
  expect_is(x[[vv$value]], "AsIs")

  # confirm a shapemap was generated and stored
  smap <- aes_map(x, "shape")
  expect_is(smap, "integer")
  expect_setequal(names(smap), x[[vv$variable]])

  # confirm shave values match manual lookup in shapemap
  expected <- I(unname(smap[x[[vv$variable]]]))
  expect_equal(x[[vv$value]], expected)
})
