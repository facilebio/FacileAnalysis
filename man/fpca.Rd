% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca-base.R, R/fpca-compare.R
\name{fpca}
\alias{fpca}
\alias{fpca.facile_frame}
\alias{compare.FacilePcaAnalysisResult}
\title{Runs a principal components analysis, the facile way.}
\usage{
fpca(
  x,
  assay_name = NULL,
  dims = 5,
  features = NULL,
  filter = "variance",
  ntop = 1000,
  row_covariates = NULL,
  col_covariates = NULL,
  batch = NULL,
  main = NULL,
  ...,
  metadata = list()
)

\method{fpca}{facile_frame}(
  x,
  assay_name = NULL,
  dims = min(5, nrow(collect(x, n = Inf)) - 1L),
  features = NULL,
  filter = "variance",
  ntop = 1000,
  row_covariates = NULL,
  col_covariates = NULL,
  batch = NULL,
  main = NULL,
  custom_key = Sys.getenv("USER"),
  ...,
  metadata = list()
)

\method{compare}{FacilePcaAnalysisResult}(x, y, run_all = TRUE, rerun = TRUE, ...)
}
\arguments{
\item{x}{a facile data container (FacileDataSet), or a \code{facile_frame}
(refer to the FacileDataStore (facile_frame) section.}

\item{assay_name}{the name of the assay to extract data from to perform the
PCA. If not specified, default assays are taken for each type of assay
container (ie. \verb{default_assay(facile container)}, \code{"counts"} for a
\code{DGEList}, \code{assayNames(SummarizedExperiment)[1L]}, etc.)}

\item{dims}{the number of PC's to calculate (minimum is 3).}

\item{features}{A feature descriptor of the features to use for the analysis.
If \code{NULL} (default), then the specified \code{filter} strategy is used.}

\item{filter}{A strategy used to identify which features to use for the
dimensionality reduction. The current (and only choice) is \code{"default"},
which takes the \code{ntop} features, sorted be decreasing variance.}

\item{ntop}{the number of features (genes) to include in the PCA. Genes are
ranked by decreasing variance across the samples in \code{x}.}

\item{row_covariates, col_covariates}{data.frames that provie meta information
for the features (rows) and samples (columns). The default is to get
these values from "the obvious places" given \code{x} (\verb{$genes} and \verb{$samples}
for a DGEList, or the sample and feature-level covariate database tables
from a FacileDataSet, for example).}

\item{batch, main}{specify the covariates to use for batch effect removal.
Refer to the \code{\link[FacileData:remove_batch_effect]{FacileData::remove_batch_effect()}} help for more information.}

\item{rerun}{when \code{rerun = TRUE} (default), the \code{fpca(x)} and \code{fpca(y)} will
be rerun over the union of the features in \code{x} and \code{y}.}
}
\value{
an fpca result
}
\description{
Performs a principal components analysis over a specified assay from the
(subset of) samples in a FacileDataStore.
}
\details{
The \code{FacilePcaAnalysisResult} produced here can be used in "the usual" ways,
ie. can be \code{viz}-ualized. \code{shine()} is 1/4th-implemented, and \code{report()}
has not been worked on yet.

Importantly / interestingly, you can shoot this result into \code{\link[=ffsea]{ffsea()}} to
perform gene set enrichment analysis over a specified dimension to identify
functional categories loaded onto differend PCs.
}
\section{Batch Correction}{

Because we assume that PCA is performed on normalized data, we leverage the
batch correction facilities provided by the \code{batch} and \code{main} parameters
in the \code{\link[FacileData:fetch_assay_data]{FacileData::fetch_assay_data()}} pipeline. If your samples have a
\code{"sex"} covariate defined, for example,  you can perform a PCA with
sex-corrected expression values like so: \code{fpca(samples, batch = "sex")}
}

\section{Features Used for PCA}{

By default, \code{fpca()} will assess the variance of all the features (genes) to
perform PCA over, and will keep the top \code{ntop} ones. This behavior is
determined by the following three parameters:
\enumerate{
\item \code{filter} determines the method by which features are selected for
analysis. Currently you can only choose \code{"variance"} (the default) or
\code{"none"}.
\item \code{features} determines the universe of features that are available for the
analysis. When \code{NULL} (default), all features for the given assay will
be loaded and filtered using the specification of the \code{filter} parameter.
If a feature descriptor is provided and \code{filter} is not specified, then
we assume that these are the exact features to run the analysis on, and
\code{filter} defaults to \code{"none"}. You may, however, intend for \code{features} to
define the universe of features to use prior to filtering, perhaps to
perform a PCA on only a certain set of genes (protein coding), but then
filter those further by variance. In this case, you will need to pass in
the feature descriptor for the universe of features you want to consider,
then \emph{explicity set \code{filter = "variance"}}.
\item \code{ntop} the default "top" number of features to take when filtering by
variance.
}
}

\section{Development Notes}{

Follow progress on implementation of \code{shine()} and \code{report()} below:
\enumerate{
\item \href{https://github.com/facilebio/FacileAnalysis/issues/12}{Implement \code{report()}}
}

Note that there are methods defined for other assay containers, like an
\code{edgeR::DGEList}, \code{limma::EList}, and \code{SummarizedExperiment}. If these are
called directly, their downstream use within the facile ecosystem isn't
yet fully supported. Development of the
\href{https://github.com/facilebio/FacileBioc}{FacileBioc package}
will address this.
}

\section{Random Things to elaborate on}{

The code here is largely inspired by DESeq2's plotPCA.

You should look at factominer:
\itemize{
\item http://factominer.free.fr/factomethods/index.html
\item http://factominer.free.fr/graphs/factoshiny.html
}
}

\section{Teaching and Tutorials}{


This looks like a useful tutorial to use when explaining the utility of
PCA analysis:
http://alexhwilliams.info/itsneuronalblog/2016/03/27/pca/

High-Dimensional Data Analysis course by Rafa Irizarry and Michael Love
https://online-learning.harvard.edu/course/data-analysis-life-sciences-4-high-dimensional-data-analysis?category[]=84&sort_by=date_added&cost[]=free
}

\section{FacileDataStore (facile_frame)}{

We enable the user to supply extra sample covariates that are not found
in the FacileDataStore associated with these samples \code{x} by adding them as
extra columns to \code{x}.

If manually provioded col_covariates have the same name as internal sample
covariates, then the manually provided ones will supersede the internals.
}

\section{Comparing PCA Results}{

We can compare two PCA results. Currently this just means we compare the
loadings of the features along each PC from fpca result \code{x} and \code{y}.
}

\examples{
efds <- FacileData::exampleFacileDataSet()

# A subset of samples ------------------------------------------------------
pca.crc <- efds |>
  FacileData::filter_samples(indication == "CRC") |>
  fpca()
if (interactive()) {
  # report(pca.crc, color_aes = "sample_type")
  shine(pca.crc)
  viz(pca.crc, color_aes = "sex")
}

# Regress "sex" out from expression data
pca.crcs <- FacileData::samples(pca.crc) |>
  fpca(batch = "sex")
if (interactive()) {
  viz(pca.crcs, color_aes = "sex")
}

# Perform PCA on only the protein coding genes
genes.pc <- features(efds) |> subset(meta == "protein_coding")
pca.crc.pc <- samples(pca.crc) |>
  fpca(features = genes.pc, filter = "variance")

pca.gdb <- pca.crc |>
  signature(dims = 1:3) |>
  result() |>
  sparrow::GeneSetDb()

# All samples --------------------------------------------------------------
pca.all <- fpca(efds)
if (interactive()) {
  viz(pca.all, color_aes = "indication", shape_aes = "sample_type")
  # report(pca.all, color_aes = "indication", shape_aes = "sample_type")
}
efds <- FacileData::exampleFacileDataSet()
p1 <- efds |>
  FacileData::filter_samples(indication == "CRC") |>
  fpca()
p2 <- efds |>
  FacileData::filter_samples(indication == "BLCA") |>
  fpca()
pcmp <- compare(p1, p2)
}
