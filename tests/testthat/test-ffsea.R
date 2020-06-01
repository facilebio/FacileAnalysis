context("Test Feature Set Enrichment Analysis (ffsea)")

if (!exists("FDS")) {
  FDS <- FacileData::exampleFacileDataSet()
}
if (!exists("gdb")) {
  gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
}

ttest.res <- FDS %>%
  FacileData::filter_samples(indication == "CRC") %>%
  flm_def(covariate = "sample_type",
          numer = "tumor", denom = "normal", batch = "sex") %>%
  fdge(method = "voom")
anova.res <- FDS %>%
  FacileData::filter_samples(indication == "CRC", sample_type == "tumor") %>%
  flm_def(covariate = "stage", batch = "sex") %>%
  fdge(method = "voom")

pca.res <- fpca(samples(anova.res))

test_that("ffsea.FacileAnalysisResult transfers all feature-level statistics to mulitGSEA result", {
  # This is really testing that the multGSEA hack currently in place that allows
  # us to pass in meta feature information and store it in the @logFC works.

  facile.gsea <- ffsea(ttest.res, gdb, methods = "cameraPR")

  # are the logFC, t-statistics, pvals, padj for each gene transferred
  # successfully?
  finfo.ttest <- ttest.res %>%
    tidy() %>%
    arrange(feature_id)

  finfo.mgres <- facile.gsea %>%
    result() %>%
    multiGSEA::logFC() %>%
    arrange(feature_id)

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
  vm <- biocbox(ttest.res)
  mgres <- multiGSEA::multiGSEA(gdb, vm, vm$design, c(-1, 1, 0),
                                "cameraPR", score.by = "logFC")

  mgres.cameraPR <- mgres %>%
    multiGSEA::result("cameraPR") %>%
    arrange(pval) %>%
    as_tibble()

  # expect_equal on tbls does not fly:
  # https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(
    select(facile.cameraPR, name, pval) %>% as.data.frame(),
    select(mgres.cameraPR, name, pval) %>% as.data.frame())
})

test_that("overrepresentation analysis (ora) works with ttest result", {
  input <- ranks(ttest.res) %>%
    tidy() %>%
    mutate(significant = padj <= 0.10)
  mgres <- multiGSEA::ora(gdb, input, selected = "significant",
                          feature.bias = "effective_length")
  mgres$name <- sub(".*;;", "", mgres$Pathway)
  facile.gsea <- ffsea(ttest.res, gdb, method = "ora",
                       biased_by = "effective_length",
                       max_padj = 0.10, min_logFC = 0)
  fres <- tidy(facile.gsea)
  expect_equal(fres$name, mgres$name)
  expect_equal(fres$pval, mgres$P.all, tolerance = 10e-4)
})

test_that("ffsea(anova_result) runs enrichment test", {
  astats <- ranks(anova.res) %>%
    tidy() %>%
    mutate(significant = padj <= 0.2)
  mgres <- multiGSEA::ora(gdb, astats, selected = "significant",
                          feature.bias = "effective_length")
  mgres$name <- sub(".*;;", "", mgres$Pathway)

  facile.gsea <- ffsea(anova.res, gdb, max_padj = 0.2,
                       biased_by = "effective_length")
  fres <- tidy(facile.gsea) %>% select(collection, name, pval, N, n, n.drawn)

  expect_equal(fres$name, mgres$name)
  expect_equal(fres$n, mgres$N)
  expect_equal(fres$pval, mgres$P.all, tolerance = 10e-4)
})

test_that("ffsea runs over dimensions of FacilePcaAnalysisResult", {
  pca1.gsea <- expect_warning({
    ffsea(pca.res, gdb, dim = 1)
  }, "Fraction.*low")
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
    "feature_id" = "feature_id", "logFC" = "score",
    "t" = "score", "score" = "score", "weight" = "weight")

  for (mname in names(check.me)) {
    fname <- check.me[mname]
    expect_equal(mgres.fstats[[mname]], pca.fstats[[fname]], info = fname)
  }
})
