if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

# High level utility ===========================================================

test_that("tidy(ANOVA) provides mean expression of all levels", {
  adef <- FDS |> 
    filter_samples(indication == "BLCA", !is.na(subtype_tcga)) |> 
    flm_def("subtype_tcga")
  adge <- fdge(adef, method = "voom")
})

test_that("tidy(simple ttest) provides mean expression of the levels tested", {
  
})

# Simpler Linear Models ========================================================

test_that("Simple fdge t-test matches explicit limma/edgeR tests", {
  mdef <- FDS |>
    filter_samples(indication == "BLCA") |>
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

  vres <- limma::lmFit(vm, vm$design) |>
    limma::eBayes() |>
    limma::topTable(coef = 2, Inf, sort.by = "none")

  expect_equal(vres$feature_id, vm_dge$feature_id)
  expect_equal(vm_dge$logFC, vres$logFC)
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
})

test_that("Simple fdge ANOVA matches explicit limma/edgeR tests", {
  mdef <- FDS |>
    filter_samples(indication == "BLCA") |>
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

  vres <- limma::lmFit(vm, vm$design) |>
    limma::eBayes() |>
    limma::topTable(coef = coefs, Inf, sort.by = "none")
  expect_equal(vm_dge$feature_id, vres$feature_id)
  expect_equal(vm_dge$pval, vres$P.Value)
  expect_equal(vm_dge$padj, vres$adj.P.Val)
  expect_equal(vm_dge$F, vres$F)
})

test_that("custom observational weights used in fdge(..., weights = W)", {
  vm.fdge <- FDS |>
    filter_samples(indication == "BLCA") |>
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal") |>
    fdge(method = "voom", with_box = TRUE) # This will calculate weights for us
  res.orig <- tidy(vm.fdge)

  vm <- biocbox(vm.fdge)

  # randomize weights and send back through
  W <- matrix(sample(vm$weights), nrow = nrow(vm))
  rownames(W) <- rownames(vm)
  colnames(W) <- colnames(vm)

  wdf <- cbind(tibble(feature_id = rownames(vm)), as.data.frame(W))

  weights <- pivot_longer(wdf, -feature_id, names_to = "sample_id",
                          values_to = "weight") |>
    separate(sample_id, "__", into = c("dataset", "sample_id"))
  
  # This is what the limma analysis looks like with custom weights
  cm <- limma::makeContrasts(tumor = tumor - normal, levels = vm$design)
  ex.res <- limma::lmFit(vm$E, vm$design, weights = W) |>
    limma::contrasts.fit(cm) |>
    limma::eBayes() |>
    limma::topTable("tumor", n = Inf, sort.by = "none")
  ex.res[["feature_id"]] <- rownames(ex.res)

  # Let's make sure user knows they can't provide custom weights for voom
  f.res <- expect_warning(
    fdge(model(vm.fdge), features = features(vm.fdge),
         method = "voom", weights = weights),
    ".*weights.*provided.*set to NULL")
  fstats <- tidy(f.res)
  expect_equal(fstats$t, res.orig$t)
  
  # We will provide our own EList with a side-laoded set of weights
  el <- vm
  el$weights <- NULL
  w.res <- fdge(model(vm.fdge), method = "limma", features = features(vm.fdge),
               biocbox = el, weights = weights)
  # we should have been warned about a different method used in the original
  # analysis vs what we are doing now.
  expect_true(any(grepl(".*method.*does not match", warnings(w.res))))
  expect_setequal(tidy(w.res)[["feature_id"]], ex.res[["feature_id"]])

  cmp <- tidy(w.res) |>
    select(feature_id, symbol, logFC, t, pval, padj) |>
    # join to stats table from manually driven weighted limma analysis
    left_join(select(ex.res, feature_id, logFC.ex = logFC, pval.ex = P.Value,
                     t.ex = t),
              by = "feature_id") |>
    # join to stats table from original voom analysis
    left_join(select(fstats, feature_id, logFC.vm = logFC, pval.vm = pval,
                     t.vm = t),
              by = "feature_id")

  expect_true(all(complete.cases(cmp))) # "full" join worked
  if (FALSE) {
    # plot similarity of logFC and tstatistics of analysis with random weignts
    # in the facile analysis vs limma-based analysis (logFC vs logFC.ex)
    plot(cmp$logFC.ex, cmp$logFC, pch = 16, col = "#3f3f3faa",
         xlab = "limma", ylab = "facile")
    abline(0, 1, col = "red")
  }
  expect_equal(cmp$logFC, cmp$logFC.ex)
  expect_equal(cmp$t, cmp$t.ex)
  expect_false(isTRUE(all.equal(cmp$t, cmp$t.vm)))
})

test_that("duplicateCorrelation is supported with voom", {
  flm <- FDS |>
    filter_samples(indication == "BLCA") |>
    flm_def(covariate = "sample_type",
            numer = "tumor",
            denom = "normal",
            block = "sex")
  expect_equal(param(flm, "block"), "sex")

  vm.fdge <- fdge(flm, method = "voom", with_box = TRUE)
  vm.box <- vm.fdge$biocbox
  assert_number(vm.box$block.corr)
  vm.res <- tidy(vm.fdge)

  # Test against two-pass voom/duplicateCorrelation mojo
  y <- samples(flm) |>
    biocbox(class = "DGEList", features = features(vm.fdge)) |>
    edgeR::calcNormFactors()
  expect_equal(nrow(y), nrow(vm.res))
  des <- model.matrix(~ 0 + sample_type, data = y$samples)

  vm <- limma::voom(y, des)
  dup <- limma::duplicateCorrelation(vm, des, block = vm$targets$sex)
  vm <- limma::voom(y, des, block = vm$targets$sex,
                    correlation = dup$consensus.correlation)
  expect_equal(vm.box$weights, vm$weights)

  dup <- limma::duplicateCorrelation(vm, des, block = vm$targets$sex)
  expect_equal(vm.box$block.corr, dup$consensus.correlation)

  ex.fit <- limma::lmFit(vm, des, block = vm$targets$sex,
                         correlation =  dup$consensus.correlation)
  ex.fit <- limma::contrasts.fit(ex.fit, c(-1, 1))
  ex.res <- ex.fit |>
    limma::eBayes() |>
    limma::topTable(n = Inf, sort.by = "none")

  expect_equal(vm.res$feature_id, rownames(ex.res))
  expect_equal(vm.res$logFC, ex.res$logFC)
  expect_equal(vm.res$t, ex.res$t)
})

test_that("a cached biocbox can be used in a call to `fdge`", {
  sf <- samples(FDS) |>
    with_sample_covariates() |>
    mutate(group = paste(indication, sample_type, sep = "_"))

  res1 <- sf |>
    flm_def("group", "BLCA_tumor", "BLCA_normal") |>
    fdge(method = "edgeR-qlf")

  expect_s4_class(biocbox(res1), "DGEList")
  
  res2 <- sf |>
    flm_def("group", "BLCA_tumor", "BLCA_normal") |>
    fdge(method = "edgeR-qlf", biocbox = biocbox(res1))
  expect_true(any(grepl("cached biobox", messages(res2))))

  t1 <- tidy(res1)
  t2 <- tidy(res2)
  expect_equal(t2, t1)
})

# Ranks and Signatures =========================================================

# Let's pre-compute and ANOVA and ttest result so we can break up our
# ranks and signature tests
TRES <- FDS |>
  filter_samples(indication == "BLCA") |>
  flm_def(covariate = "sample_type",
                 numer = "tumor",
                 denom = "normal") |>
  fdge(method = "voom")
ARES <- FDS |>
  filter_samples(indication == "BLCA") |>
  flm_def(covariate = "stage", batch = "sex") |>
  fdge(method = "voom")

test_that("ttest ranks and signatures generated correctly", {
  # Signed Ranks
  eranks.signed <- TRES |>
    tidy() |>
    arrange(desc(logFC))
  ranks.signed <- TRES |>
    ranks(signed = TRUE, rank_by = "logFC") |>
    tidy()
  expect_equal(ranks.signed[["feature_id"]], eranks.signed[["feature_id"]])

  # Unsigned Ranks
  eranks.unsigned <- TRES |>
    tidy() |>
    arrange(pval)
  ranks.unsigend <- TRES |>
    ranks(signed = FALSE, rank_by = "pval") |>
    tidy()
  expect_equal(ranks.unsigend[["feature_id"]], eranks.unsigned[["feature_id"]])
})

test_that("ranks returns anova features in expected order", {
  eranks <- ARES |>
    tidy() |>
    arrange(desc(.data[["F"]]))
  anova.ranks <- ARES |>
    ranks() |>
    tidy()
  expect_equal(anova.ranks[["feature_id"]], eranks[["feature_id"]])

  # signed anova ranks fires a warning but returns a ranking
  r <- expect_warning(ranks(ARES, signed = TRUE), "only.*unsigned")
  expect_equal(tidy(r)[["feature_id"]], eranks[["feature_id"]])
})

test_that("t-test signatures generated correcty", {
  # we'll compare these to parts of the correctly-generated ranks, which were
  # tested above.
  sig.signed <- signature(TRES, signed = TRUE) |> tidy()
  sig.unsigned <- signature(TRES, signed = FALSE) |> tidy()

  nup.signed <- sum(sig.signed$direction == "up")
  expect_equal(nup.signed, nrow(sig.signed) / 2)
  expect_equal(nup.signed, sum(sig.signed$direction == "down"))

  joint <- intersect(sig.unsigned$feature_id, sig.signed$feature_id)
  expect_gt(length(joint), nup.signed - 1)
})


