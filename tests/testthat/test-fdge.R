context("Differential Gene Expression and GSEA")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Simple fdge t-test matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "edgeR-qlf")
  qlf_dge <- result(qlf_test)

  bbox <- biocbox(qlf_test)
  y <- result(bbox)
  expect_equal(nrow(y), nrow(qlf_dge))
  expect_equal(rownames(y), qlf_dge$feature_id)

  design <- model.matrix(~ sample_type, y$samples)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  qres <- edgeR::glmQLFTest(qfit, coef = 2)
  qres <- edgeR::topTags(qres, n = Inf, sort.by = "none")
  qres <- edgeR::as.data.frame.TopTags(qres)

  expect_equal(qres$feature_id, qlf_dge$feature_id)
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$logFC, qres$logFC)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom")
  vm_dge <- result(vm_test)

  vm <- limma::voom(y, y$design)
  vm.facile <- biocbox(vm_test)$biocbox
  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = 2, Inf, sort.by = "none")

  expect_equal(vres$feature_id, vm_dge$feature_id)
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$logFC, vres$logFC)
})

test_that("Simple fdge ANOVA matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "stage", fixed = "sex")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "edgeR-qlf")
  qlf_dge <- result(qlf_test)

  bbox <- biocbox(qlf_test)
  y <- result(bbox)
  design <- model.matrix(~ stage + sex, y$samples)
  coefs <- grep("stage", colnames(design))

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  qres <- edgeR::glmQLFTest(qfit, coef = coefs)
  qres <- edgeR::topTags(qres, n = Inf)
  qres <- edgeR::as.data.frame.TopTags(qres)
  qres <- qres[qlf_dge$feature_id,]
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$F, qres$F)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom")
  vm_dge <- result(vm_test)

  vm <- limma::voom(y, y$design)
  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = coefs, Inf)
  vres <- vres[vm_dge$feature_id,]
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$F, vres$F)
})

# Let's pre-compute and ANOVA and ttest result so we can break up our
# ranks and signature tests
TRES <- FDS %>%
  filter_samples(indication == "BLCA") %>%
  fdge_model_def(covariate = "sample_type",
                 numer = "tumor",
                 denom = "normal") %>%
  fdge(method = "voom")
ARES <- FDS %>%
  filter_samples(indication == "BLCA") %>%
  fdge_model_def(covariate = "stage", fixed = "sex") %>%
  fdge(method = "voom")

test_that("ttest ranks and signatures generated correctly", {
  # Signed Ranks
  eranks.signed <- TRES %>%
    result() %>%
    arrange(desc(logFC))
  ranks.signed <- TRES %>%
    ranks(signed = TRUE) %>%
    result()
  expect_equal(ranks.signed[["feature_id"]], eranks.signed[["feature_id"]])

  # Unsigned Ranks
  eranks.unsigned <- TRES %>%
    result() %>%
    arrange(pval)
  ranks.unsigend <- TRES %>%
    ranks(signed = FALSE) %>%
    result()
  expect_equal(ranks.unsigend[["feature_id"]], eranks.unsigned[["feature_id"]])
})

test_that("ranks returns anova features in expected order", {
  eranks <- ARES %>%
    result() %>%
    arrange(pval)
  anova.ranks <- ARES %>%
    ranks() %>%
    result()
  expect_equal(anova.ranks[["feature_id"]], eranks[["feature_id"]])

  # signed anova ranks fires a warning but returns a ranking
  r <- expect_warning(ranks(ARES, signed = TRUE), "only.*unsigned")
  expect_equal(result(r)[["feature_id"]], eranks[["feature_id"]])
})

test_that("t-test signatures generated correcty", {
  # we'll compare these to parts of the correctly-generated ranks, which were
  # tested above.
  sig.signed <- signature(TRES, signed = TRUE)

  sig.unsinged <- signature(TRES, signed = FALSE)
})


