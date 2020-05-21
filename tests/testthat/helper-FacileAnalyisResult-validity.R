# Validaty of analysis result objects ------------------------------------------
# These functions were developed to test validity of FacileAnalysisResult
# objects after fsave/fload cycle

test_valid_fdge <- function(x, fds, ...) {
  # note: currently not using the reference `fds` for anything
  expect_class(fds, "FacileDataStore")

  expect_is(x, "FacileDgeAnalysisResult")
  expect_s3_class(fds(x), "FacileDataStore")
  expect_s3_class(fds(model(x)), "FacileDataStore")
  expect_s3_class(fds(samples(x)), "FacileDataStore")
}

test_valid_ffsea <- function(x, fds, ...) {
  # note: currently not using the reference `fds` for anything
  expect_class(fds, "FacileDataStore")

  expect_is(x, "FacileFseaAnalysisResult")
  mgres <- result(x)
  expect_is(mgres, "MultiGSEAResult")

  parent <- param(x, "x")
  if (is(parent, "FacileDgeAnalysisResult")) {
    test_valid_fdge(parent, fds)
  }
}
