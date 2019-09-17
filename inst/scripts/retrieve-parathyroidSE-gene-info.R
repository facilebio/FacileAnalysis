library(SummarizedExperiment)

data("parathyroidGenesSE", package = "parathyroidSE")

bm <- loadNamespace("biomaRt")
mart <- bm$useMart(
  host = "feb2014.archive.ensembl.org",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl")
mart.info <- bm$getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = rownames(parathyroidGenesSE),
  mart = mart)

write.csv(mart.info, "inst/extdata/parathyroidSE-gene-info.csv",
          row.names=FALSE)
write.csv(mart.info, "inst/extdata/airway-gene-info.csv",
          row.names=FALSE)
