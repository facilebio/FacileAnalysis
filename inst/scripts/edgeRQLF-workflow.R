# Helper functions for the edgeRQLF workflow vignette

#' Encapsulates DGEList creation code from the RnaSeqGeneEdgeRQL workflow
#'
#' As the assembly of a DGEList isn't really something we care to teach,
#' Who wants to put all of this code in a vignette?
assemble_edgeRQLF_DGEList <- function(tmpdir = tempdir()) {
  # ............................................................... Targets File
  targets <- system.file("extdata", "targets.txt",
                         package = "RnaSeqGeneEdgeRQL")
  targets <- read.delim(targets, stringsAsFactors=FALSE)
  targets$group <- factor(paste(targets$CellType, targets$Status, sep="."))

  # ..................................................................... counts
  counts.fn <- file.path(tmpdir, "GSE60450_Lactation-GenewiseCounts.txt.gz")
  if (!file.exists(counts.fn)) {
    FileURL <- paste(
      "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
      "format=file",
      "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
      sep="&")
    download.file(FileURL, method = "libcurl", counts.fn)
  }

  GenewiseCounts <- read.delim(counts.fn, row.names = "EntrezGeneID")
  colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)

  # ...................................................................... genes
  a <- loadNamespace("AnnotationDbi")
  genes <- GenewiseCounts[, 1L, drop=FALSE]
  genes[["symbol"]] <- a$mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                rownames(GenewiseCounts),
                                keytype = "ENTREZID", column = "SYMBOL")
  edgeR::DGEList(GenewiseCounts[, -1L], genes = genes,
                 samples = transform(targets, GEO = NULL, SRA = NULL))
}
