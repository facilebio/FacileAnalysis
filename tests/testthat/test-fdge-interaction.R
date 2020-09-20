context("Differential Gene Expression: interaction / compare")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

anova.all <- samples(FDS) %>%
  with_sample_covariates(c("indication", "sample_type", "sex")) %>%
  mutate(group = paste(indication, sample_type, sep = "_")) %>%
  flm_def("group", batch = "sex") %>%
  fdge()
tvn.blca <- samples(anova.all) %>%
  flm_def("group", "BLCA_tumor", "BLCA_normal", batch = "sex") %>%
  fdge(features = features(anova.all))
tvn.crc <- samples(anova.all) %>%
  flm_def("group", "CRC_tumor", "CRC_normal", batch = "sex") %>%
  fdge(features = features(anova.all))
tvn.compare <- compare(tvn.blca, tvn.crc)
istats <- tidy(tvn.compare)

test_that("ttest compare() over same covariate/numer/denom, disjoint samples", {
  y.all <- anova.all %>%
    samples() %>%
    biocbox("DGEList", features = features(anova.all))
  y.all <- suppressWarnings(edgeR::calcNormFactors(y.all))

  des <- model.matrix(~ 0 + group + sex, y.all$samples)
  colnames(des) <- sub("group", "", colnames(des))

  ctr <- limma::makeContrasts(
    blca = BLCA_tumor - BLCA_normal,
    crc = CRC_tumor - CRC_normal,
    interaction = (BLCA_tumor - BLCA_normal) - (CRC_tumor - CRC_normal),
    levels = des)
  vm <- limma::voom(y.all, des)
  fit <- limma::lmFit(vm, vm$des) %>%
    limma::contrasts.fit(ctr) %>%
    limma::eBayes()

  lres <- fit %>%
    limma::topTable("interaction", n = Inf, sort.by = "none") %>%
    select(feature_id, logFC, t, pval = P.Value, padj = adj.P.Val)
  expect_set_equal(istats$feature_id, lres$feature_id)
  lres <- lres[istats$feature_id,]

  expect_equal(istats$logFC, lres$logFC)
  expect_equal(istats$t, lres$t)
  expect_equal(istats$pval, lres$pval)
})

test_that("interaction logFC is correct when full tests are not run", {
  approx <- compare(tvn.blca, tvn.crc, .run_interaction = FALSE)
  astats <- tidy(approx)

  expect_equal(istats$feature_id, astats$feature_id)
  expect_equal(istats$logFC.x, astats$logFC.x)
  expect_equal(istats$logFC.y, astats$logFC.y)
  expect_equal(istats$pval.x, astats$pval.x)
  expect_equal(istats$pval.y, astats$pval.y)

  expect_equal(istats$logFC, astats$logFC)
})

if (FALSE) {
test_that("ttest compare() over same covariate, different numer/denom", {
  # subset to cancer type, (stage I vs stage II) vs (stage III vs stage IV)
  crc.samples <- filter_samples(FDS, indication == "CRC")
  stage.2v1 <- crc.samples %>%
    flm_def("stage", "II", "I") %>%
    fdge()
  stage.4v3 <- crc.samples %>%
    flm_def("stage", "IV", "III") %>%
    fdge()
  cmp.stage <- compare(stage.2v1, stage.4v3)
})
}
