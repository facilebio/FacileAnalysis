#' Save and load FacileAnalysisResult objects to/from disk.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Saving and loading FacileAnalysisResults can be tricky because they rely
#' on a chain of parent objects that created them, all of which have a reference
#' to the FacileDataStore that they came from.
#'
#' Currenly, the FacileDataStore is removed from every internal object when
#' before it is serialized via `fsave()`. The FacileDataStore is then restored
#' to "just where it needs to be" (using strong assumptions of where it needs
#' to be, instead of logging where it was removed from in `fsave()`) when
#' brought back to life via `fload()`.
#'
#' @details
#' fsave creates a directory with two elements, the actually analysis result
#' (a qs2 object), and a `meta.yaml` that has metadata related to the 
#' FacileAnalysisResult object being saved. The details in the `meta.yaml`
#' file will include the `metadata()` from `x`, and details about what
#' FacileDataStore was used to serialize the result.
#' 
#' metadata about the analysis result will be stored under the `result:`
#' section, and the `datastore:` section will hold the elemnts that were
#' previously associated with the `fsave_info` attribute from the v1
#' *.qs object
#' 
#' OLD:
#' An `fsave_info` attribute will be attached to the serialized object before
#' saving which will be a list that holds information required to materialize
#' the analysis result successfully. Minimally this will include info about how
#' to reconstitute and attach the FacileDataStore that the object was generated
#' from.
#'
#' Elements of `fsave_info` list so far include:
#'
#' - `fds_class`: A string indicating the class of the FacileDataStore that the
#'    result `x` was generated from.
#' - `fds_path`: The path on the filesystem that teh FacileDataStore can be
#'   found. Currently we expect this to be a FacileDataSet
#' - `fds_anno_dir`: The annotation directory for the FacileDataSet
#' - `fds_add`: A descriptor that indicate which elements in the serialized
#'   object should have the fds attached to them, besides itself.
#'
#' @export
#' @rdname serialize
#'
#' @param x A `FacileAnalysisResult` object.
#' @param file A directory named `file` will be created to serialize the bits
#'   required to reload this analysis result in another session. If `file`
#'   already exists, this function will throw an error.
#' @param with_fds Serialize the FacileDataStore with the object? Default is
#'   `FALSE` and you should have a good reason to change this behavior.
#' @return `NULL` for `fsave`, the (correctly sublcassed) `FacileAnalysisResult`
#'   for `fload`, along with `FacileDataStore` (`fds`) attached
fsave <- function(x, file, with_fds = FALSE, ...) {
  UseMethod("fsave", x)
}

#' @noRd
#' @export
fsave.FacileDgeAnalysisResult <- function(x, file, with_fds = FALSE,
                                          with_biocbox = FALSE, ...) {
  if (!with_biocbox) x[["biocbox"]] <- NULL
  NextMethod()
}

#' @noRd
#' @export
fsave.FacileAnalysisResult <- function(x, file, with_fds = FALSE, ...) {
  lifecycle::signal_stage("experimental", "fsave()")
  fds. <- assert_class(fds(x), "FacileDataStore")

  checkmate::assert_directory_exists(dirname(file), "w")
  if (file.exists(file)) stop("Destination `", file, "` already exists")
  if (!dir.create(file)) stop("Failed to create output directory: ", file)
  
  x.info <- metadata(x)
  x.info$class <- class(x)

  fds.info <- list(class = class(fds.))
  if ("FacileDataSet" %in% fds.info[["class"]]) {
    fds.info[["path"]] <- fds.[["parent.dir"]]
    fds.info[["anno_dir"]] <- fds.[["anno.dir"]]
  }
  
  info <- list(
    result = x.info,
    datastore = fds.info
  )
  yaml::write_yaml(info, file.path(file, "meta.yaml"))
  
  x <- unfds(x)
  qs2::qs_save(x, file.path(file, "result.qs2"))
  invisible(file)
}


#' @rdname serialize
#' @export
#' @param fds The `FacileDataStore` the object was run on.
fload <- function(x, fds = NULL, anno = NULL, with_fds = TRUE, ...) {
  checkmate::assert_directory_exists(x, "r")
  res.fn <- checkmate::assert_file_exists(file.path(x, "result.qs2"))
  meta.fn <- checkmate::assert_file_exists(file.path(x, "meta.yaml"))
  meta <- yaml::read_yaml(meta.fn)
  res <- qs2::qs_read(res.fn)
  assert_class(res, "FacileAnalysisResult")

  fds.info <- meta[["datastore"]]
  if (is.null(fds.info)) {
    stop("Meta information about connected FacileDataStore not found, ",
         "did you save this object using the FacileAnalysis::fsave() function?")
  }
  assert_list(fds.info, names = "unique")
  fds.class <- checkmate::assert_character(fds.info[["class"]], min.len = 1L)

  if (with_fds) {
    if (is.null(fds) || checkmate::test_string(fds)) {
      # Implicit FacileDataStores only work with FacileDataSet for now, since
      # it needs to exist in a serialized form on the filesystem anyway, so
      # there is somewhere we can retrieve it from
      if (!"FacileDataSet" %in% fds.class) {
        stop("FacileDataSet required if `fds` not provided")
      }
      if (is.null(fds)) {
        fds <- fds.info[["path"]]
      }
      if (!test_directory_exists(fds, "r")) {
        stop("No path to linked FacileDataStore object found in `file`,
           please pass in a value for `fds` explicitly")
      }
      if (is.null(anno)) {
        anno <- fds.info[["anno_dir"]]
      }
      fds <- FacileData::FacileDataSet(fds, anno.dir = anno)
    }
    if (!is(fds, "FacileDataStore")) {
      stop("We expected a FacileDataStore by now")
    }
    # common.class <- intersect(fds.info[["fds_class"]], class(fds))
    # if (length(common.class) == 0L) {
    #   stop("It doesn't look like the FacileDataStore you are trying to ",
    #        "associate with this result is the one that was used")
    # }
    
    xs <- samples(res)
    fs <- dplyr::collect(samples(fds), n = Inf)
    xmissed <- dplyr::anti_join(xs, fs, by = c("dataset", "sample_id"))
    if (nrow(xmissed) > 0L) {
      stop(nrow(xmissed), " samples missing from the provided faciledatastore")
    }
    
    res <- refds(res, fds)
  }

  res
}
