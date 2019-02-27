context("Differential Gene Expression and GSEA")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("Simple fdge t-test matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "sample_type",
                   numer = "tumor",
                   denom = "normal")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "qlf", gsea = NULL)
  qlf_dge <- qlf_test$dge

  y <- fdge_biocbox(mdef, assay_name = "rnaseq", method = "qlf")
  design <- model.matrix(~ sample_type, y$samples)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  qres <- edgeR::glmQLFTest(qfit, coef = 2)
  qres <- edgeR::topTags(qres, n = Inf)
  qres <- edgeR::as.data.frame.TopTags(qres)
  qres <- qres[qlf_test$dge$feature_id,]
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$logFC, qres$logFC)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom", gsea = NULL)
  vm_dge <- vm_test$dge

  vm <- limma::voom(y, y$design)
  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = 2, Inf)
  vres <- vres[vm_dge$feature_id,]
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$logFC, vres$logFC)
})

test_that("Simple fdge ANOVA matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    fdge_model_def(covariate = "stage", fixed = "sex")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "qlf", gsea = NULL)
  qlf_dge <- qlf_test$dge

  y <- fdge_biocbox(mdef, assay_name = "rnaseq", method = "qlf")
  design <- model.matrix(~ stage + sex, y$samples)
  coefs <- grep("stage", colnames(design))

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  qres <- edgeR::glmQLFTest(qfit, coef = coefs)
  qres <- edgeR::topTags(qres, n = Inf)
  qres <- edgeR::as.data.frame.TopTags(qres)
  qres <- qres[qlf_test$dge$feature_id,]
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$F, qres$F)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom", gsea = NULL)
  vm_dge <- vm_test$dge

  vm <- limma::voom(y, y$design)
  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = coefs, Inf)
  vres <- vres[vm_dge$feature_id,]
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$F, vres$F)
})
