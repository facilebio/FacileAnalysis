library("testthat")
library("checkmate")
library("dplyr")
library("FacileData")
library("FacileAnalysis")

FDS <- FacileData::exampleFacileDataSet()

test_check("FacileAnalysis")

