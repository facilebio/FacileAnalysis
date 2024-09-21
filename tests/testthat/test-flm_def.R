if (!exists("afds")) afds <- FacileData::an_fds()
if (!exists("xsamples")) {
  xsamples <- afds |> 
    FacileData::filter_samples(cell_abbrev %in% c("CNT", "EC"))
}

test_that("flm_def errors when numerator and denominator are the same", {
  expect_error({
    flm_def(xsamples, "cell_abbrev", "CNT", "CNT")
  }, "numer.*denom.*same")
})