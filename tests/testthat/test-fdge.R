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
  qlf_dge <- dge(qlf_test)

  y <- biocbox(qlf_test)
  expect_equal(nrow(y), nrow(qlf_dge))
  expect_equal(rownames(y), qlf_dge$feature_id)

  design <- model.matrix(~ sample_type, y$samples)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qres <- edgeR::glmQLFit(y, y$design, robust = TRUE) %>%
    edgeR::glmQLFTest(coef = 2) %>%
    edgeR::topTags(n = Inf, sort.by = "none") %>%
    edgeR::as.data.frame.TopTags()

  if (FALSE) {
    plot(qlf_dge$logFC, qres$logFC, pch = 16, col = "#3f3f3f33")
    abline(0, 1, col = "red")
  }

  expect_equal(qres$feature_id, qlf_dge$feature_id)
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$logFC, qres$logFC)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom", gsea = NULL)
  vm_dge <- dge(vm_test)

  vm <- biocbox(vm_test)
  expect_equal(nrow(vm), nrow(vm_dge))
  expect_equal(rownames(vm), vm_dge$feature_id)

  design <- model.matrix(~ sample_type, vm$targets)

  vres <- limma::lmFit(vm, design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = 2, Inf, sort.by = "none")
  expect_equal(vres$feature_id, vm_dge$feature_id)

  if (FALSE) {
    plot(vm_dge$logFC, vres$logFC, pch = 16, col = "#3f3f3f33")
    abline(0, 1, col = "red")
  }

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
  qlf_dge <- dge(qlf_test)

  y <- biocbox(qlf_test)
  expect_equal(nrow(y), nrow(qlf_dge))
  expect_equal(rownames(y), qlf_dge$feature_id)

  design <- model.matrix(~ stage + sex, y$samples)
  coefs <- grep("stage", colnames(design))

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  qres <- edgeR::glmQLFit(y, y$design, robust = TRUE) %>%
    edgeR::glmQLFTest(coef = coefs) %>%
    edgeR::topTags(n = Inf, sort.by = "none") %>%
    edgeR::as.data.frame.TopTags()

  expect_equal(qres$feature_id, qlf_dge$feature_id)
  expect_equal(qlf_dge$pval, qres$PValue)
  expect_equal(qlf_dge$padj, qres$FDR)
  expect_equal(qlf_dge$F, qres$F)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom", gsea = NULL)
  vm_dge <- dge(vm_test)

  vm <- biocbox(vm_test)
  expect_equal(nrow(vm), nrow(vm_dge))
  expect_equal(rownames(vm), vm_dge$feature_id)

  design <- model.matrix(~ stage + sex, vm$targets)
  coefs <- grep("stage", colnames(design))

  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = coefs, Inf, sort.by = "none")
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$F, vres$F)
})
