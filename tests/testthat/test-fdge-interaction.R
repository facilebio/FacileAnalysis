context("Differential Gene Expression: interaction / compare")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

test_that("ttest compare() over same covariate/numer/denom, disjoint samples", {
  tvn.blca <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def("sample_type", "tumor", "normal", batch = "sex") %>%
    fdge()
  tvn.crc <- FDS %>%
    filter_samples(indication == "CRC") %>%
    flm_def("sample_type", "tumor", "normal", batch = "sex") %>%
    fdge()
  tvn.cmp <- expect_warning(compare(tvn.blca, tvn.crc), "ALPHA")

  all.samples <- samples(tvn.crc) %>%
    bind_rows(samples(tvn.blca)) %>%
    as_facile_frame(FDS) %>%
    transmute(dataset, sample_id,
              group = paste(dataset, sample_type, sep = "_"))

  y.all <- as.DGEList(all.samples, feature_ids = features(tvn.cmp))
  y.all <- suppressWarnings(calcNormFactors(y.all)) # partial match of `p` to `probs`
  des <- model.matrix(~ 0 + group + sex, y.all$samples)
  colnames(des) <- sub("group", "", colnames(des))
  ctr <- makeContrasts(
    testme = (BLCA_tumor - BLCA_normal) - (COAD_tumor - COAD_normal),
    levels = des)
  vm <- limma::voom(y.all, des)
  ires <- limma::lmFit(vm, vm$des) %>%
    limma::contrasts.fit(ctr) %>%
    limma::eBayes() %>%
    limma::topTable("testme", n = Inf) %>%
    select(feature_id, logFC, t, pval = P.Value, padj = adj.P.Val)

  expect_set_equal(tidy(tvn.cmp)$feature_id, ires$feature_id)
  cmp <- tidy(tvn.cmp) %>%
    select(feature_id, symbol, logFC, t, pval, padj) %>%
    inner_join(ires, by = "feature_id")
  expect_equal(cmp$logFC.x, cmp$logFC.y)
  expect_equal(cmp$t.x, cmp$t.y)
  expect_equal(cmp$pval.x, cmp$pval.y)
})

test_that("ttest compare() over same covariate, different numer/denom", {
  # subset to cancer type, (stage I vs stage II) vs (stage III vs stage IV)
  crc.samples <- filter_samples(FDS, indication == "CRC")
  stage.2v1 <- crc.samples %>%
    flm_def("stage", "II", "I") %>%
    fdge()
  stage.4v3 <- crc.samples %>%
    flm_def("stage", "IV", "III") %>%
    fdge()
  cmp.stage <- expect_warning(compare(stage.2v1, stage.4v3), "ALPHA")
})
