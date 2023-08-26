if (!exists("FDS")) {
  FDS <- FacileData::exampleFacileDataSet()
}
if (!exists("gdb")) {
  gdb <- sparrow::getMSigGeneSetDb("h", "human", id.type = "entrez")
}

ttest.res <- FDS |>
  FacileData::filter_samples(indication == "CRC") |>
  flm_def(covariate = "sample_type",
          numer = "tumor", denom = "normal", batch = "sex") |>
  fdge(method = "voom")

gsea.res <- ffsea(ttest.res, gdb, "cameraPR")


test_that("unfds/refds removes and restores 'facility' of fdge results", {
  fds. <- fds(ttest.res)
  expect_s3_class(fds., "FacileDataStore")

  # strip
  uttest <- unfds(ttest.res)
  expect_null(fds(uttest))

  # reconstitute
  rttest <- refds(uttest, fds.)
  test_valid_fdge(rttest, fds.)
})

test_that("unfds/refds removes and restores 'facility' of ffsea results", {
  fds. <- fds(ttest.res)
  expect_s3_class(fds., "FacileDataStore")

  # strip
  ugsea <- unfds(gsea.res)
  expect_null(fds(ugsea))
  expect_null(param(ugsea, "x") |> fds())
  expect_warning(fds(samples(ugsea)), "no.*datastore.*found", ignore.case = TRUE)

  # reconstitute
  rgsea <- refds(ugsea, fds.)
  test_valid_ffsea(rgsea, fds.)
})

test_that("fsave/fload brings back the attached FacileDataStore", {
  fds. <- fds(gsea.res)
  expect_s3_class(fds., "FacileDataStore")

  withr::with_tempfile("fc", fileext = ".rds", {
    fsave(gsea.res, fc)
    test_valid_ffsea(fload(fc), fds.)
    test_valid_ffsea(fload(fc, fds = fds.), fds.)
    # no datastore requested
    # res.nofds <- fload(fn, fds = fds., with_fds = FALSE)
  })
})
