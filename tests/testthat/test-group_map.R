if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()

# More tidy like analyses
test_that("group_map fdge is identical to lapply/sample subset", {
  exs <- samples(FDS) |> 
    with_sample_covariates() |> 
    filter(indication == "BLCA")
  
  dge.list <- sapply(c("m", "f"), function(sx) {
    exs |> 
      filter(.data$sex == sx) |> 
      flm_def("sample_type", "tumor", "normal") |> 
      fdge(method = "voom", metadata = list(sex = sx))
  }, simplify = FALSE)
  
  dge.map <- exs |> 
    group_by(sex) |> 
    group_map.(~ {
      .x |> 
        flm_def("sample_type", "tumor", "normal") |> 
        fdge(method = "voom", metadata = list(sex = .y$sex))
    })
  
  for (sex in names(dge.list)) {
    cmp <- full_join(
      tidy(dge.list[[sex]]),
      tidy(dge.map[[sex]]), 
      by = "feature_id", suffix = c(".list", ".map"))
    expect_equal(
      cmp$logFC.list, 
      cmp$logFC.map, 
      info = paste0("logFC do not match for `", sex, "`"))
    expect_equal(
      cmp$pval.list, 
      cmp$pval.map,
      info = paste0("pvalues do not match for `", sex, "`"))
  }

})