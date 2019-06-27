
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileAnalysis

<!-- badges: start -->

[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/facileverse/FacileAnalysis.svg?branch=master)](https://travis-ci.org/facileverse/FacileAnalysis)
[![Codecov test
coverage](https://codecov.io/gh/facileverse/FacileAnalysis/branch/master/graph/badge.svg)](https://codecov.io/gh/facileverse/FacileAnalysis?branch=master)
<!-- badges: end -->

The FacileAnalysis package defines a set of analysis tasks over genomic
data in a modular fashion, which requires these data to be stored in a
container that implements the FacileData API (aka a `FacileDataStore`).

The over-arching goal of developing analyses within this framework is to
enable quick and effortless interrogation of genomic data by enabling a
**hybrid** interactive and code-driven approach to data analysis.
Analyses can be either completely code-driven, GUI driven, or some mix
of the two.

To achieve this goal, analysis modules break down a general analysis
task into smaller constituent steps, the results of which can:

1.  be piped together (`%>%`) to perform a complete analysis;
2.  be explored at **different levels of interactivity**, via the
    `shine()`, `viz()`, and `report()` methods; and
3.  act as starting (or reference) points for **the next** analysis
    task, via the `compare()`, `ranks()`, and `signature()` methods.

You can find a more in-depth overview and philosophy of this package in
the [FacileAnalysis vignette](FacileAnalysis.html).

The analyses implemented within this package are listed below, with
links to vignettes that describe their functionality in more detail:

  - `fpca`: [Interactive Principal Components
    Analyses](FacileAnalysis-fpca.html)
  - `fdge`: [Interactive Differential Gene Expression
    Analysis](FacileAnalysis-fdge.html)
  - `fsea`: [Interactive (Gene) Set Enrichment
    Analysis](FacileAnalysis-fsea.html)

## Example Analysis: Differential Expression

We’ll include an exemplar differential expression analysis here in order
to touch on some of the guiding principles of an analysis modules
mentioned above.

We have first defined a complete differential expression analysis
(`fdge`), by breaking it down into the following steps:

1.  identifying the subset of samples from a `FacileDataStore` to
    perform the analysis over;
2.  defining the model (design) over the samples, ie. the covariate to
    test and which covariates to used for fixed/batch effects
3.  defining the parameters to run the statistical test, ie.
    1)  the assay from the samples to run the analysis on
    2)  the differential expression pipeline to use (limma, limma-trend,
        voom ,etc)
    3)  advanced options, like the threshold to test against
        (limma/treat), or whether to incorporate sample level weights

The example below identifies genes differentially expression between
tumor and normal samples (and controlling for sex) in the “BLCA”
indication of the example TCGA dataset included in the
[FacileData](https://facileverse.github.io/FacileData) package.

``` r
library(FacileData)
library(FacileAnalysis)
efds <- exampleFacileDataSet()

# Step 1: define the samples implicated in our test
samples <- filter_samples(efds, indication == "BLCA")

# Step 2: define the model (design) for the test
model <- fdge_model_def(samples, covariate = "sample_type",
                        numer = "tumor", denom = "normal",
                        fixed = "sex")

# Step 3: configure the options to run the test, which include the assay that
#         holds the data used to test and the statistical method/framework we
#         should use to perform the test
vdge <- fdge(model, assay_name = "rnaseq", method = "voom")
```

Perhaps you prefer the
[edgeR/QLF](https://f1000research.com/articles/5-1438) analysis
framework, instead? No problem, we only need to tweak one of the
parameters in the last step of the pipeline:

``` r
qdge <- fdge(model, assay_name = "rnaseq", method = "edgeR-qlf")
```

… or DESeq2, perhaps? No problem, [we accept pull
requests](https://github.com/denalitherapeutics/FacileAnalysis/pulls)
:-)

**Note** that `fdge` is running over the RNA-seq data stored in the
“rnaseq” assay. If the datastore had a `"CNV"` assay matrix, one could
find regions of differential copy number between tumor and normal like
so:

``` r
cnv.dge <- fdge(model, assay_name = "CNV", method = "limma-trend")
```

We are using a `limma-trend` pipeline there, which might not be
appropriate. Developers can develop their own `method` for differential
copy number analysis and plug it in there.

## Interrogation of Analyses at Different Levels of Interactivity

There are a number of S3 methods a `FacileAnalysisResult` needs to
define in order to be complete. The `shine()`, `viz()`, and `report()`
methods allow the end-user to interrogate (and report) the results of an
analysis at different levels of interactivity.

Let’s take a look at how these work over a `FacileDgeAnalysisResult`

### Deep Interaction via shine()

The `shine(aresult, ...)` method provides the richest interactive view
over a `FacileAnalysisResult` by launching a [shiny
gadget](https://shiny.rstudio.com/articles/gadgets.html) that enables
the end-user to fully interrogate the results.

``` r
shine(vdge)
```

<img src="man/figures/fdge-shine.png" width="50%" />

### Interactive Graphics via viz()

The `viz(aresult, ...)` methods leverage
[htmlwidgets](https://www.htmlwidgets.org) to create a
JavaScript-powered interactive view of the analysis result, which is
detached from a running R-process.

``` r
viz(vdge)
```

### Reporting Results via report()

While the output of the `viz()` functions can be used directly in
Rmarkdown reports, `report(aresult, ...)` is meant to create a “more
complete” (perhaps multi-htmlwidget) view over the result that can be
more suitable for inclusion into an Rmarkdown report.

``` r
cap <- paste(
  "Differential expression of tumor vs normal samples, in the TCGA bladder"
  "cancer indication (BLCA)")
report(vdge, caption = cap)
```

## Hybrid Analyses

The same differential expression analysis that created the `vdge` object
above can be performed entirely interactively, with the same results.

We can either start the analysis from the same predefined set of
samples, but define the linear model and testing framework to use
interactively, by launching a the **f**acile **d**ifferential **g**ene
**e**xpression gadget:

``` r
vdge2 <- fdgeGadget(samples)
```

Or we can perform the whole thing via a GUI which lets us select the
subset of samples and run the differential expression analysis without
using any code at all:

``` r
vdge3 <- fdgeGadget(efds)
```

Assuming the same filtering and testing strategies were selected in the
GUI using the `fdgeGadget` calls above, the objects they return will all
be equivalent to the `vdge` result, which was entirely generated
programmatically.

## Installation

``` r
# These packages are only available via github for now
# install.packages("devtools")
devtools::install_github("facileverse/FacileAnalysis")
```

## Resources

  - [facileverse.github.io/FacileAnalysis](https://facileverse.github.io/FacileAnalysis)
    (online documentation and vignettes)
  - [Ask a
    question](https://groups.google.com/forum/#!forum/facileverse)
    (Google Group Q\&A Site)
  - [Open an
    issue](https://github.com/facileverse/FacileAnalysis/issues) (GitHub
    issues for bug reports and feature requests)
