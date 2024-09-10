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
#' @param file A path to a file to save. If filename ends with '*.rds', then
#'   object will be serialized with `saveRDS`. If filename ends with '*.qs',
#'   then `qs::qsave` will be used. In teh future, we should support
#'   `base::connections`.
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
  ext <- tolower(tools::file_ext(file))
  if (ext == "rds") {
    save.fn <- saveRDS
  } else if (ext == "qs") {
    reqpkg("qs")
    save.fn <- qs::qsave
  } else {
    stop("Unknown filetype extension to save: ", ext)
  }

  lifecycle::signal_stage("experimental", "fsave()")
  fds. <- assert_class(fds(x), "FacileDataStore")

  fds.info <- list(fds_class = class(fds.))
  if ("FacileDataSet" %in% fds.info[["fds_class"]]) {
    fds.info[["fds_path"]] <- fds.[["parent.dir"]]
    fds.info[["fds_anno_dir"]] <- fds.[["anno.dir"]]
  }

  x <- unfds(x)
  attr(x, "fsave_info") <- fds.info
  save.fn(x, file)
  invisible(file)
}


#' @rdname serialize
#' @export
#' @param fds The `FacileDataStore` the object was run on.
fload <- function(x, fds = NULL, anno = NULL, with_fds = TRUE, ...) {
  if (test_string(x)) {
    ext <- tolower(tools::file_ext(x))
    if (ext == "rds") {
      read.fn <- readRDS
    } else if (ext == "qs") {
      reqpkg("qs")
      read.fn <- qs::qread
    } else {
      stop("Unknown filetype extension to load: ", ext)
    }
    x <- read.fn(x)
  }
  assert_class(x, "FacileAnalysisResult")

  fds.info <- attr(x, "fsave_info")
  if (is.null(fds.info)) {
    stop("Meta information about connected FacileDataStore not found, ",
         "did you save this object using the FacileAnalysis::fsave() function?")
  }
  assert_list(fds.info, names = "unique")
  assert_subset(c("fds_class"), names(fds.info))
  fds.class <- assert_character(fds.info[["fds_class"]], min.len = 1L)

  if (with_fds) {
    if (is.null(fds) || test_string(fds)) {
      # Implicit FacileDataStores only work with FacileDataSet for now, since
      # it needs to exist in a serialized form on the filesystem anyway, so
      # there is somewhere we can retrieve it from
      if (!"FacileDataSet" %in% fds.class) {
        stop("FacileDataSet required if `fds` not provided")
      }
      if (is.null(fds)) {
        fds <- fds.info[["fds_path"]]
      }
      if (!test_directory_exists(fds, "r")) {
        stop("No path to linked FacileDataStore object found in `file`,
           please pass in a value for `fds` explicitly")
      }
      if (is.null(anno)) {
        anno <- fds.info[["fds_anno_dir"]]
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
    
    xs <- samples(x)
    fs <- dplyr::collect(samples(fds), n = Inf)
    xmissed <- dplyr::anti_join(xs, fs, by = c("dataset", "sample_id"))
    if (nrow(xmissed) > 0L) {
      stop(nrow(xmissed), " samples missing from the provided faciledatastore")
    }
    
    x <- refds(x, fds)
  }

  x
}
