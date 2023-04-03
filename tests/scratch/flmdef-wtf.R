devtools::install_local("~/facilebio/FacileShine", force = TRUE)

library(FacileData)
devtools::load_all(".")

fds <- FacileDataSet("/datasets/datastores/FacileDenaliCynoDataSet")

flm <- fds |>
  filter_samples(study_name == "DST-215") |>
  flm_def("group", "AL258_PSEG__mgkg60", "nonTV_TREM2__mgkg50")
res <- fdge(flm, filter_universe = features(fds) |> filter(meta == "protein_coding"))

r2 <- fdgeGadget(samples(flm))
