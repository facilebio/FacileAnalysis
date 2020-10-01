context("Interaction models: dge and gsea")

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

# Differential Expression ------------------------------------------------------
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

# GSEA -------------------------------------------------------------------------
if (!exists("gdb.h")) {
  gdb.h <- multiGSEA::getMSigGeneSetDb("h", "human", id.type = "entrez")
}

if (!exists("gdb.go")) {
  if (!exists("gdb.go.bp")) {
    gdb.go.bp <- multiGSEA::getMSigGeneSetDb("c5", "human", id.type = "entrez")
    gdb.go.bp <- gdb.go.bp[geneSets(gdb.go.bp)$collection == "GO_BP"]
  }
  set.seed(0xBEEF)
  .pselected <- 0.05
  .gk <- sample(c(TRUE, FALSE), length(gdb.go.bp), replace = TRUE,
                prob = c(.pselected, 1 - .pselected))
  .some.gs <- c(
    "DNA_STRAND_ELONGATION",
    "DNA_STRAND_ELONGATION_INVOLVED_IN_DNA_REPLICATION",
    "GENE_SILENCING_BY_RNA",
    "MULTICELLULAR_ORGANISMAL_SIGNALING",
    "REGULATION_OF_CELL_DIFFERENTIATION",
    "REGULATION_OF_OSTEOBLAST_DIFFERENTIATION",
    "FLAVONOID_GLUCURONIDATION",
    "HOMOPHILIC_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",
    "URONIC_ACID_METABOLIC_PROCESS",
    "XENOBIOTIC_GLUCURONIDATION")

  .gk <- .gk | geneSets(gdb.go.bp)$name %in% .some.gs
  gdb.go <- gdb.go.bp[.gk]
}

if (FALSE) {
  cviz <- viz(tvn.compare)
}

test_that("ffsea on interaction stats works like a normal ttest ffsea", {
  max.padj <- 0.10
  min.logFC <- 1
  res.gsea.h <- ffsea(tvn.compare, gdb.h, methods = c("cameraPR", "ora"),
                      # tweak a default parameter
                      rank_by = "t",
                      # these are default values for params, but being explicit
                      max_padj = max.padj, min_logFC = min.logFC,
                      group_by = "direction")
  exp.gsea.h <- tidy(tvn.compare) %>%
    mutate(direction = ifelse(logFC > 0, "up", "down"),
           significant = abs(logFC) >= min.logFC & padj <= max.padj) %>%
    ffsea(gdb = param(res.gsea.h, "gdb"),
          methods = param(res.gsea.h, "methods"),
          rank_by = param(res.gsea.h, "rank_by"),
          rank_order = "descending",
          select_by = "significant",
          group_by = param(res.gsea.h, "group_by"))
  res.names <- c("cameraPR", "ora", "ora.up", "ora.down")
  for (rname in res.names) {
    res.stats <- tidy(res.gsea.h, rname)
    exp.stats <- tidy(exp.gsea.h, rname)
    expect_equal(res.stats, exp.stats, info = rname)
  }
})

test_that("quadrant analysis runs ora over each quadrant", {
  quad.res <- ffsea(tvn.compare, gdb.go, type = "quadrants")
  expect_setequal(names(quad.res$quadrants), c("both", "x", "y"))
})
