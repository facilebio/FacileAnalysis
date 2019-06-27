context("Test (Gene) Set Enrichment Analysis (fsea)")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("fsea,cameraPR works like multiGSEA call", {
  gdb <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")

  # GSEA the facile way
  facile.ttest <- FDS %>%
    FacileData::filter_samples(indication == "CRC") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor", denom = "normal", fixed = "sex") %>%
    fdge(method = "voom")
  facile.gsea <- fsea(facile.ttest, gdb, method = "cameraPR")
  facile.mgres <- result(facile.gsea)
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
