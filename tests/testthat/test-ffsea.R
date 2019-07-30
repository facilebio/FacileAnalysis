context("Test Feature Set Enrichment Analysis (ffsea)")

if (!exists("FDS")) {
  FDS <- FacileData::exampleFacileDataSet()
}
if (!exists("gdb")) {
  gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
}

ttest.res <- FDS %>%
  FacileData::filter_samples(indication == "CRC") %>%
  fdge_model_def(covariate = "sample_type",
                 numer = "tumor", denom = "normal", fixed = "sex") %>%
  fdge(method = "voom")
anova.res <- FDS %>%
  FacileData::filter_samples(indication == "CRC", sample_type == "tumor") %>%
  fdge_model_def(covariate = "stage", fixed = "sex") %>%
  fdge(method = "voom")

pca.res <- fpca(samples(anova.res))

test_that("ffsea.FacileAnalysisResult transfers all feature-level statistics to mulitGSEA result", {
  # This is really testing that the multGSEA hack currently in place that allows
  # us to pass in meta feature information and store it in the @logFC works.

  facile.gsea <- ffsea(ttest.res, gdb, method = "cameraPR")

  # are the logFC, t-statistics, pvals, padj for each gene transferred
  # successfully?
  finfo.ttest <- ttest.res %>%
    tidy() %>%
    arrange(feature_id)

  finfo.mgres <- facile.gsea %>%
    result() %>%
    multiGSEA::logFC() %>%
    arrange(featureId) %>%
    rename(feature_id = "featureId")

  expect_equal(nrow(finfo.mgres), nrow(finfo.ttest))
  test.cols <- c("feature_id", "t", "logFC", "pval", "padj", "CI.L", "CI.R")
  for (cname in test.cols) {
    expect_equal(finfo.mgres[[cname]], finfo.ttest[[cname]], info = cname)
  }
})

test_that("cameraPR call through ffsea works like multiGSEA call", {
  # GSEA the facile way
  facile.gsea <- ffsea(ttest.res, gdb, method = "cameraPR")

  facile.cameraPR <- facile.gsea %>%
    tidy(name = "cameraPR") %>%
    arrange(pval)

  # GSEA the "traditional" multiGSEA way
  bbox <- biocbox(ttest.res)
  vm <- result(bbox)
  mgres <- multiGSEA::multiGSEA(gdb, vm, vm$design, c(-1, 1, 0),
                                "cameraPR", score.by = "logFC")

  mgres.cameraPR <- mgres %>%
    multiGSEA::result("cameraPR") %>%
    arrange(pval) %>%
    as.tbl()

  # expect_equal on tbls does not fly:
  # https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(
    select(facile.cameraPR, name, pval) %>% as.data.frame(),
    select(mgres.cameraPR, name, pval) %>% as.data.frame())
})

test_that("ffsea runs over dimensions of FacilePcaAnalysisResult", {
  pca1.gsea <- ffsea(pca.res, gdb, dim = 1)
  mgres <- result(pca1.gsea)

  # check that pc-stuff are in features of gsea result.
  # currently (multiGSEA_v0.12.8), the mgres@logFC $logFC and $t columns will
  # be loaded with the "score" of the fpca feature rankings
  pca.fstats <- ranks(pca.res, signed = TRUE, dims = 1) %>% tidy()
  mgres.fstats <- mgres %>%
    multiGSEA::logFC() %>%
    arrange(desc(score))
  expect_equal(nrow(mgres.fstats), nrow(pca.fstats))

  # feature_id, logFC, t, and score should all be the same
  # name is multGSEA column, value is ranks column
  check.me <- c(
    "featureId" = "feature_id", "logFC" = "score",
    "t" = "score", "score" = "score", "weight" = "weight")

  for (mname in names(check.me)) {
    fname <- check.me[mname]
    expect_equal(mgres.fstats[[mname]], pca.fstats[[fname]], info = fname)
  }
})

if (FALSE) {
  # are there really no hallmark pathways differentially expressed between
  # tumor and normal samples?
  y <- filter_samples(FDS, indication == "CRC") %>% as.DGEList()
  des <- model.matrix(~ sample_type + sex, y$samples)
  y <- y[filterByExpr(y, des),]
  vm <- voom(y, des)
  mg2 <- multiGSEA(gdb, vm, vm$design, "sample_typetumor",
                   c("camera", "cameraPR"))
  # guess not
}
