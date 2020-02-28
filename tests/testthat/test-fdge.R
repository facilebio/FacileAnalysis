context("Differential Gene Expression")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

# Simpler Linear Models ========================================================

test_that("Simple fdge t-test matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "edgeR-qlf",
                   with_box = TRUE)
  qlf_dge <- tidy(qlf_test)
  bb <- qlf_test$biocbox

  y <- biocbox(qlf_test)

  expect_is(y, "DGEList")
  expect_equal(nrow(y), nrow(bb))
  expect_equal(rownames(y), rownames(bb))
  expect_equal(colnames(y), colnames(bb))
  expect_equal(y$counts, bb$counts, check.attributes = FALSE)
  expect_equal(y$samples$lib.size, bb$samples$lib.size)
  expect_equal(y$samples$norm.factors, bb$samples$norm.factors)

  design <- model.matrix(~ sample_type, y$samples)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)

  expect_equal(y$trend.method, bb$trend.method)
  expect_equal(y$common.dispersion, bb$common.dispersion)

  qfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)


  qres <- edgeR::glmQLFTest(qfit, coef = 2)
  qres <- edgeR::topTags(qres, n = Inf, sort.by = "none")
  qres <- edgeR::as.data.frame.TopTags(qres)

  expect_equal(qres$feature_id, qlf_dge$feature_id)
  expect_equal(qres$PValue, qlf_dge$pval)
  expect_equal(qres$FDR, qlf_dge$padj)
  expect_equal(qres$logFC, qlf_dge$logFC)

  # Test voom
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom")
  vm_dge <- result(vm_test)
  expect_set_equal(vm_dge[["feature_id"]], qlf_dge[["feature_id"]])

  vm <- limma::voom(y, y$design)
  vm.facile <- biocbox(vm_test)
  expect_equal(vm.facile$weights, vm$weights)

  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = 2, Inf, sort.by = "none")

  expect_equal(vres$feature_id, vm_dge$feature_id)
  expect_equal(vm_dge$logFC, vres$logFC)
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
})

test_that("Simple fdge ANOVA matches explicit limma/edgeR tests", {
  mdef <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def(covariate = "stage", batch = "sex")

  # Test edgeR quasilikelihood
  qlf_test <- fdge(mdef, assay_name = "rnaseq", method = "edgeR-qlf")
  qlf_dge <- tidy(qlf_test)

  y <- biocbox(qlf_test)
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
  vm_test <- fdge(mdef, assay_name = "rnaseq", method = "voom", with_box = TRUE)
  vm_dge <- tidy(vm_test)
  vmf <- vm_test$biocbox

  vm <- limma::voom(y, y$design)
  expect_equal(vm$weights, vmf$weights)

  vres <- limma::lmFit(vm, vm$design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = coefs, Inf, sort.by = "none")
  expect_equal(vm_dge$feature_id, vres$feature_id)
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$F, vres$F)
})

test_that("custom observational weights used in fdge(..., weights = W)", {
  # TODO: test weighted linear models!
  vm.fdge <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal") %>%
    fdge(method = "voom", with_box = TRUE) # This will calculate weights for us

  vm <- biocbox(vm.fdge)

  # randomize weights and send back through
  W <- matrix(sample(vm$weights), nrow = nrow(vm))
  rownames(W) <- rownames(vm)
  colnames(W) <- colnames(vm)

  wdf <- cbind(tibble(feature_id = rownames(vm)), as.data.frame(W))

  weights <- pivot_longer(wdf, -feature_id, names_to = "sample_id",
                          values_to = "weight") %>%
    separate(sample_id, "__", into = c("dataset", "sample_id"))

  cm <- limma::makeContrasts(tumor = tumor - normal, levels = vm$design)
  ex.res <- limma::lmFit(vm$E, vm$design, weights = W) %>%
    limma::contrasts.fit(cm) %>%
    limma::eBayes() %>%
    limma::topTable("tumor", n = Inf, sort.by = "none") %>%
    mutate(feature_id = rownames(.))

  f.res <- expect_warning(
    fdge(model(vm.fdge), features = features(vm.fdge),
         method = "voom", weights = weights),
    ".*weights.*replaced")
  fstats <- tidy(f.res)

  expect_setequal(tidy(vm.fdge)[["feature_id"]], ex.res[["feature_id"]])
  expect_setequal(fstats[["feature_id"]], ex.res[["feature_id"]])

  cmp <- tidy(vm.fdge) %>%
    select(feature_id, symbol, pval, padj) %>%
    left_join(select(ex.res, feature_id, logFC.limma = logFC, pval.limma = P.Value),
              by = "feature_id") %>%
    left_join(select(fstats, feature_id, logFC.f2 = logFC, pval.f2 = pval),
              by = "feature_id")

  expect_true(all(complete.cases(cmp))) # "full" join worked
  expect_false(isTRUE(all.equal(cmp[["logFC"]], cmp[["logFC.limma"]])))
  expect_equal(cmp[["logFC.f2"]], cmp[["logFC.limma"]])
  expect_equal(cmp[["pval.f2"]], cmp[["pval.limma"]])
})

# Ranks and Signatures =========================================================

# Let's pre-compute and ANOVA and ttest result so we can break up our
# ranks and signature tests
TRES <- FDS %>%
  filter_samples(indication == "BLCA") %>%
  flm_def(covariate = "sample_type",
                 numer = "tumor",
                 denom = "normal") %>%
  fdge(method = "voom")
ARES <- FDS %>%
  filter_samples(indication == "BLCA") %>%
  flm_def(covariate = "stage", batch = "sex") %>%
  fdge(method = "voom")

test_that("ttest ranks and signatures generated correctly", {
  # Signed Ranks
  eranks.signed <- TRES %>%
    tidy() %>%
    arrange(desc(logFC))
  ranks.signed <- TRES %>%
    ranks(signed = TRUE) %>%
    tidy()
  expect_equal(ranks.signed[["feature_id"]], eranks.signed[["feature_id"]])

  # Unsigned Ranks
  eranks.unsigned <- TRES %>%
    tidy() %>%
    arrange(pval)
  ranks.unsigend <- TRES %>%
    ranks(signed = FALSE) %>%
    tidy()
  expect_equal(ranks.unsigend[["feature_id"]], eranks.unsigned[["feature_id"]])
})

test_that("ranks returns anova features in expected order", {
  eranks <- ARES %>%
    tidy() %>%
    arrange(pval)
  anova.ranks <- ARES %>%
    ranks() %>%
    tidy()
  expect_equal(anova.ranks[["feature_id"]], eranks[["feature_id"]])

  # signed anova ranks fires a warning but returns a ranking
  r <- expect_warning(ranks(ARES, signed = TRUE), "only.*unsigned")
  expect_equal(tidy(r)[["feature_id"]], eranks[["feature_id"]])
})

test_that("t-test signatures generated correcty", {
  # we'll compare these to parts of the correctly-generated ranks, which were
  # tested above.
  sig.signed <- signature(TRES, signed = TRUE)

  sig.unsinged <- signature(TRES, signed = FALSE)
})


