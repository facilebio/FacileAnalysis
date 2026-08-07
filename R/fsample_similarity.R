#' Run correlation analysis across a set of samples.
#'
#' This is a typical QC step for RNA-seq data, for instance, where you want
#' to show how well samples are correlated with each other, hopefully to see
#' that the samples from withint the same experimental group are better
#' correlated with each other than samples from other gruops
#'
#' @export
#' @param x a facile_frame of samples
#' @param method "pearson" (default), or "spearman"
#' @param covariate if a string, then pull off the column from `x`
#' @param features what features to use?
#' @param filter default to "variance"
#' @param ntop default to top 1000
#' @examples
#' afds <- FacileData::an_fds()
#' xs.all <- FacileData::samples(afds) |> FacileData::with_sample_covariates()
#' xs <- xs.all |> dplyr::filter(cond == "DKD")
#' ss <- fsample_similarity(xs)
#' 
#' ctypes <- unique(xs$cell_abbrev)
#' cell.cols <- scales::pal_brewer(palette = "Set3")(length(ctypes))
#' names(cell.cols) <- ctypes
#' 
#' viz(ss, "corr", annotate = "cell_abbrev", color_map = cell.cols)
#' viz(ss, "dist", annotate = "cell_abbrev", color_map = cell.cols)
fsample_similarity <- function(
    x,
    assay_name = NULL,
    cor_method = "pearson",
    dist_method = "euclidean",
    use = "everything",
    features = NULL,
    filter = "variance",
    ntop = 1000,
    col_covariates = NULL,
    batch = NULL,
    main = NULL,
    ...,
    metadata = list()
) {
  UseMethod("fsample_similarity", x)
}

#' @noRd
#' @export
fsample_similarity.facile_frame <- function(
    x,
    assay_name = NULL,
    cor_method = "pearson",
    dist_method = "euclidean",
    use = "everything",
    features = NULL,
    filter = "variance",
    ntop = 1000,
    col_covariates = NULL,
    batch = NULL,
    main = NULL,
    ...,
    metadata = list()
) {
  cor_method <- match.arg(cor_method, .cor.methods())
  dist_method <- match.arg(dist_method, .dist.methods())
  use <- match.arg(use, .na.methods())
  
  # this is taken from FacileAnalysis::fpca.facile_frame -----------------------
  .fds <- assert_class(FacileData::fds(x), "FacileDataStore")
  FacileData::assert_sample_subset(x)
  x <- FacileData::collect(x, n = Inf)
  
  messages <- character()
  warnings <- character()
  errors <- character()
  
  if (is.null(assay_name)) {
    assay_name <- FacileData::default_assay(.fds)
  }
  
  # This should do the batch effect removal
  if (FacileViz::unselected(batch)) {
    batch <- NULL
  }
  if (FacileViz::unselected(main)) {
    main <- NULL
  }
  
  # Subset the samples from `x` that have values for the given assay under test.
  xs <- FacileData::filter_by_assay_support(x, assay_name)
  dropped <- FacileData::samples(xs, dropped = TRUE)
  ndropped <- nrow(dropped)
  if (ndropped > 0) {
    msg <- paste(
      ndropped,
      "samples have no",
      assay_name,
      "data. These samples will ",
      "be removed for downstream analysis."
    )
    warnings <- c(warnings, msg)
    x <- xs
    
    if (nrow(x) == 0L) {
      stop("No samples left after filtering against assay availability")
    }
  }
  
  dat <- FacileData::biocbox(
    x,
    class = "list",
    assay_name = assay_name,
    features = features,
    sample_covariates = col_covariates,
    normalized = TRUE,
    log = TRUE,
    # batch = batch, main = main,
    ...,
    metadata = metadata
  )

  if (!is.null(features) && missing(filter)) {
    filter <- if (!missing(ntop)) "variance" else "none"
  }
  
  # end fpca code --------------------------------------------------------------
  
  out <- fsample_similarity(
    dat[["assay_data"]],
    cor_method = cor_method,
    dist_method = dist_method,
    use = use,
    features = features,
    filter = filter,
    ntop = ntop,
    col_covariates = dat[["samples"]],
    batch = batch,
    main = main,
    ...,
    metadata = metadata
  )
  
  out$samples <- x
  out$features <- dplyr::as_tibble(dat[["features"]][out$features,,drop=FALSE])
  out
}

#' @noRd
#' @export
fsample_similarity.matrix <- function(
    x,
    cor_method = "pearson",
    dist_method = "euclidean",
    use = "everything",
    features = NULL,
    filter = "variance",
    ntop = 1000,
    col_covariates = NULL,
    batch = NULL,
    main = NULL,
    ...,
    metadata = list()
) {
  cor_method <- match.arg(cor_method, .cor.methods())
  dist_method <- match.arg(dist_method, .dist.methods())
  use <- match.arg(use, .na.methods())
  
  messages <- character()
  warnings <- character()
  errors <- character()
  
  if (!is.null(features)) {
    features <- extract_feature_id(features)
    features <- unique(features)
  }
  if (!is.null(features) && missing(filter)) {
    filter <- "none"
  }
  assert_choice(filter, c("variance", "none"))
  
  if (is.null(rownames(x))) {
    stop("No rownames() on assay matrix: this is unusual")
    rownames(x) <- as.character(seq(nrow(x)))
  }
  
  if (is(col_covariates, "data.frame")) {
    assert_true(ncol(x) == nrow(col_covariates))
    assert_character(rownames(col_covariates))
    assert_true(all(colnames(x) == rownames(col_covariates)))
  }
  
  if (FacileViz::unselected(batch)) {
    batch <- NULL
  }
  if (FacileViz::unselected(main)) {
    main <- NULL
  }
  
  # Some assays let NA's sneak in, if we observe any here let's just drop
  # the analyte and "pray for the best". We will introduce an NA fill policy
  # soon.
  isna <- which(is.na(x), arr.ind = TRUE)
  if (nrow(isna) > 0L) {
    rm.na <- unique(isna[, 1L])
    warning(
      "Removing ",
      length(rm.na),
      " / ",
      nrow(x),
      " features due to NA values",
      immediate. = TRUE
    )
    if (is.character(features)) {
      features <- setdiff(features, rownames(x)[rm.na])
    }
    x <- x[-rm.na, , drop = FALSE]
  }
  
  if (!is.null(batch)) {
    x <- remove_batch_effect(x, col_covariates, batch = batch, main = main, ...)
  }
  
  if (!is.null(features)) {
    take <- match(features, rownames(x))
    if (any(is.na(take))) {
      stop(
        "The fcormap (fpca) filtering strategy only allows you to specify rownames ",
        "(features) to use for PCA, and the ones you specified do not ",
        "exactly match rownames(x), see:\n  ",
        "https://github.com/facilebio/FacileAnalysis/issues/20"
      )
    }
  }
  
  if (filter == "none") {
    take <- seq(nrow(x))
    ntop <- nrow(x)
  } else if (filter == "variance") {
    rv <- matrixStats::rowVars(x, na.rm = TRUE)
    take <- head(order(rv, decreasing = TRUE), ntop)
    ntop <- length(take)
  }
  
  xx <- x[take, , drop = FALSE]
  
  distances <- dist(t(xx), method = dist_method)
  cormat <- cor(xx, method = cor_method, use = use)
  out <- list(
    correlation = cormat,
    distance = distances,
    features = rownames(xx),
    metadata = metadata,
    params = list(
      cor_method = cor_method,
      dist_method = dist_method,
      use = use,
      filter = filter,
      ntop = ntop,
      batch = batch,
      main = main
    )
  )
  class(out) <- c("FacileSimilarityResult", "FacileAnalysisResult")
  out
}

#' @noRd
#' @export
samples.FacileSimilarityResult <- function(x, ...) {
  x$samples
}

#' @noRd
#' @importFrom FacileAnalysis result
#' @export
result.FacileSimilarityResult <- function(
    x,
    name = c("correlation", "distance"),
    ...
) {
  name <- match.arg(name)
  x[[name]]
}

#' @noRd
#' @export
viz.FacileSimilarityResult <- function(
    x,
    name = c("distance", "correlation"),
    annotate = NULL,
    corrplot_method = "circle",
    hclust_method = NULL,
    color_map = NULL,
    strip_dataset_prefix = TRUE,
    corr_coef_size = 8, # set to 0 if you don't want to include
    title = NULL,
    legend_position = c("bottom", "side"),
    ...
) {
  name <- match.arg(name)
  legend_position <- match.arg(legend_position)
  xx <- result(x, name)
  xm <- if (name == "distance") as.matrix(xx) else xx
  if (is.null(title)) title <- name
  checkmate::assert_string(title)

  if (is.null(hclust_method)) {
    hclust_method <- if (name == "distance") "ward.D2" else "complete"
  }
  hclust_method <- match.arg(hclust_method, c("ward.D2", "complete", "none"))
  
  adf <- NULL
  if (checkmate::test_character(annotate)) {
    xs <- samples(x)
    amissing <- setdiff(annotate, colnames(xs))
    xs <- FacileData::with_sample_covariates(xs, amissing)
    amissing <- setdiff(annotate, colnames(xs))
    if (length(amissing)) {
      stop(
        "Could not find covariates to annotate viz: ",
        paste(amissing, collapse = ",")
      )
    }
    adf <- as.data.frame(xs[, annotate, drop = FALSE])
    rownames(adf) <- paste(xs$dataset, xs$sample_id, sep = "__")
    
    stopifnot(setequal(rownames(adf), rownames(xm)))
    adf <- adf[rownames(xm), , drop = FALSE]
  }

  if (strip_dataset_prefix) {
    rownames(xm) <- sub(".*?__", "", rownames(xm))
    colnames(xm) <- rownames(xm)
    rownames(adf) <- rownames(xm)
  }
  
  if (checkmate::test_character(color_map) && length(annotate) == 1L) {
    color_map <- list(colors = color_map)
    names(color_map) <- annotate
  }
  if (!checkmate::test_list(color_map, names = "unique")) {
    warning("Can't parse color_map correctly, setting to NULL")
    color_map <- NULL
  }
  
  col.anno <- ComplexHeatmap::HeatmapAnnotation(
    df = adf,
    col = color_map
  )
  
  if (legend_position == "bottom") {
    hm_legend_param <- list(legend_direction = "horizontal")
  } else {
    hm_legend_param <- list(legend_direction = "verticalks")
  }
  
  if (name == "distance") {
    dcols <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(200)
    # out <- ComplexHeatmap::pheatmap(
    #   xm,
    #   clustering_distance_rows = xx,
    #   clustering_distance_cols = xx,
    #   col = dcols,
    #   annotation_col = adf,
    #   annotation_colors = color_map
    # )
    
    out <- ComplexHeatmap::Heatmap(
      xm,
      col = circlize::colorRamp2(
        # c(0, quantile(as.vector(xm), 0.98)),
        c(0, max(xm)),
        c("#08306B", "#F7FBFF")
      ),
      name = "Distance",
      clustering_distance_rows = xx,
      clustering_distance_columns = xx,
      top_annotation = col.anno,
      column_title = title,
      heatmap_legend_param = hm_legend_param
    )
  } else {
      # cor.pal <- colorRampPalette(rev(cor.cols))(200)
      # out <- corrplot::corrplot(
      #   xm,
      #   col = cor.pal,
      #   method = corrplot_method,
      #   hclust.method = "ward.D2",
      #   tl.col = "grey30"
      # )

    out <- ComplexHeatmap::Heatmap(
      xm,
      col = circlize::colorRamp2(
        c(-1, 0, 1),
        c("navy", "white", "firebrick")
      ),
      name = "Correlation",
      rect_gp = grid::gpar(type = "none"), # Turn off default background rectangles

      cell_fun = function(j, i, x, y, width, height, fill) {
        # Get absolute value for diameter scaling
        r <- abs(xm[i, j])
        # Draw background white grid cells for clean layout
        grid::grid.rect(
          x,
          y,
          width,
          height,
          gp = grid::gpar(col = "#eeeeee", fill = "white")
        )
        # Draw circle with diameter proportional to the absolute correlation
        # Max diameter is min(width, height) when abs(r) = 1
        grid::grid.circle(
          x, y, 
          # r = r * 0.5 * min(grid::unit.c(width, height)),
          # from corrplot
          #             circles = asp_rescale_factor * 0.9 * abs(DAT) ^ 0.5 / 2,
          r = r^0.5/2 * min(grid::unit.c(width, height)),
          gp = grid::gpar(fill = fill, col = NA)
        )
        # include correlation coefficient
        if (checkmate::test_int(corr_coef_size, lower = 1)) {
          grid::grid.text(
            sprintf("%.2f", r),
            x,
            y,
            gp = grid::gpar(fontsize = corr_coef_size)
          )
        }
      },
      
      cluster_rows = hclust_method != "none",
      clustering_method_rows = hclust_method,

      cluster_columns = hclust_method != "none",
      clustering_method_columns = hclust_method,

      top_annotation = col.anno,
      column_title = title,
      heatmap_legend_param = hm_legend_param
    )
  }
  
  out
}

.dist.methods <- function() {
  c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
}

.cor.methods <- function() {
  c("pearson", "kendall", "spearman")
}

.na.methods <- function() {
  c(
    "all.obs",
    "complete.obs",
    "pairwise.complete.obs",
    "everything",
    "na.or.complete"
  )
}
