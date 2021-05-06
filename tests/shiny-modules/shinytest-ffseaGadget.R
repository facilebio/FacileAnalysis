library(FacileData)

devtools::load_all(".")
options(facile.log.level.fshine = "trace")
options(facile.log.level.fanalysis = "trace")

efds <- FacileData::exampleFacileDataSet()
gdb <- sparrow::getMSigGeneSetDb("h", "human", "entrez")

pca.crc <- efds %>%
  FacileData::filter_samples(indication == "CRC") %>%
  fpca()
pca.gsea <- ffseaGadget(pca.crc, gdb)

if (FALSE) {
  pca.gsea <- ffsea(pca.crc, gdb, dim = 1, methods = c("cameraPR", "fgsea"))
  shine(pca.gsea)
}


ttest.res <- efds %>%
  FacileData::filter_samples(indication == "CRC") %>%
  flm_def(covariate = "sample_type",
                 numer = "tumor", denom = "normal", batch = "sex") %>%
  fdge(method = "voom")
ttest.gsea <- ffseaGadget(ttest.res, gdb)

stage.anova <- efds %>%
  FacileData::filter_samples(indication == "BLCA") %>%
  flm_def(covariate = "stage", batch = "sex") %>%
  fdge(method = "voom")
anova.gsea <- ffseaGadget(stage.anova, gdb)

