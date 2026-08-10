volcano_test_data <- function() {
  data.frame(
    symbol = c("A", "Dup", "Dup", NA_character_, ""),
    feature_id = paste0("gene_", letters[1:5]),
    logFC = c(-2, -1, 0, 1, 2),
    pval = c(0.001, 0.01, 0.5, 0.02, 1),
    FDR = c(0.005, 0.04, 0.6, 0.08, 1),
    stringsAsFactors = FALSE
  )
}

test_that("volcano metric helpers fall back and transform y values", {
  y <- c("p-value" = "pval", "FDR" = "FDR")
  dat <- volcano_test_data()

  expect_equal(.fvolcano_y_metric(NULL, y), "pval")
  expect_equal(.fvolcano_y_metric("bad", y), "pval")
  expect_equal(.fvolcano_y_metric("FDR", y), "FDR")
  expect_equal(.fvolcano_y_label("FDR"), "-log10(FDR)")
  expect_equal(.fvolcano_y_values(dat, "pval"), -log10(dat$pval))
})

test_that("volcano labels are disambiguated and limited", {
  dat <- volcano_test_data()
  choices <- fvolcanoLabelChoices(
    dat,
    c("symbol", "feature_id"),
    "feature_id"
  )
  expect_equal(
    names(choices),
    c("A", "Dup (gene_b)", "Dup (gene_c)", "gene_d", "gene_e")
  )

  y.values <- .fvolcano_y_values(dat, "pval")
  labels <- .fvolcano_label_data(
    dat,
    y.values,
    c("gene_a", "gene_b", "gene_c"),
    "logFC",
    "feature_id",
    c("symbol", "feature_id"),
    label_limit = 2L
  )
  expect_equal(labels$feature_id, c("gene_a", "gene_b"))
})

test_that("volcano point styles classify cutoff and selected points", {
  dat <- volcano_test_data()
  y.values <- .fvolcano_y_values(dat, "pval")
  style <- .fvolcano_point_style(
    dat,
    y.values,
    x.cutoff = 1,
    y.cutoff = 1.5,
    selected = "gene_c",
    x = "logFC",
    key = "feature_id",
    label_limit = 20L
  )

  expect_equal(style$color[1], .fvolcano_point_palette$left)
  expect_equal(style$color[3], .fvolcano_point_palette$selected)
  expect_equal(style$line.width[c(1, 3)], c(0.5, 0.5))
  expect_equal(style$line.width[5], 0)
})

test_that("volcano axis ranges handle empty and constant values", {
  empty <- data.frame(logFC = numeric(), yaxis = numeric())
  expect_equal(.fvolcano_axis_ranges(empty, "logFC")$x, c(-0.58, 0.68))

  constant <- data.frame(logFC = c(1, 1), yaxis = c(2, 2))
  ranges <- .fvolcano_axis_ranges(constant, "logFC")
  expect_lt(ranges$x[1], 1)
  expect_gt(ranges$x[2], 1)
  expect_lt(ranges$y[1], 2)
  expect_gt(ranges$y[2], 2)
})

test_that("fvolcanoPlot builds a plotly widget", {
  dat <- volcano_test_data()
  plt <- fvolcanoPlot(
    dat,
    y = "pval",
    labels = "gene_a",
    x_cutoff = 1,
    y_cutoff = 1.5,
    hover = NULL,
    event_source = "volcano_test"
  )

  expect_s3_class(plt, "plotly")
  annotations <- unlist(
    lapply(plt$x$layoutAttrs, function(x) length(x$annotations))
  )
  expect_true(any(annotations == 1L))
})
