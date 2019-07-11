context("Test Feature Set Enrichment Analysis (ffsea)")

if (!exists("FDS")) {
  FDS <- FacileData::exampleFacileDataSet()
}
if (!exists("gdb")) {
  gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
}

test_that("cameraPR call through ffsea works like multiGSEA call", {
  # GSEA the facile way
  facile.ttest <- FDS %>%
    FacileData::filter_samples(indication == "CRC") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor", denom = "normal", fixed = "sex") %>%
    fdge(method = "voom")
  facile.gsea <- ffsea(facile.ttest, gdb, method = "cameraPR")
  facile.mgres <- result(facile.gsea, "object")
  facile.cameraPR <- multiGSEA::result(facile.mgres, "cameraPR")

  # GSEA the "traditional" multiGSEA way
  bbox <- biocbox(facile.ttest)
  vm <- result(bbox)
  mgres <- multiGSEA::multiGSEA(gdb, vm, vm$design, c(-1, 1, 0),
                                "cameraPR", score.by = "logFC")

  # First test that the differential expression statistics are the same, so
  # that any failure on the following gsea test is confined to the GSEA and
  # not having tested the wrong dge question
  expect_equal(result(facile.ttest)$logFC, multiGSEA::logFC(mgres)$logFC)

  mgres.cameraPR <- multiGSEA::result(mgres, "cameraPR")
  expect_equal(
    select(facile.cameraPR, name, pval),
    select(mgres.cameraPR, name, pval))
})

test_that("ffsea runs over dimensions of FacilePcaAnalysisResult", {
  pca.crc <- FacileData::exampleFacileDataSet() %>%
    FacileData::filter_samples(indication == "CRC") %>%
    fpca()
  pca1.gsea <- ffsea(pca.crc, gdb, pc = 1)
  mgres <- result(pca1.gsea)
  # TODO: test something here ffsea,pca
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
