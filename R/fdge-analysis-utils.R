#' Tally the number of up/down regulated genes across a series of
#' comparisons.
#'
#' @seealso [plot_dge_summary()]
#' @export
#' @param comps a data.frame of comparisons that were run, with a `$label`
#'   column
#' @param dge.res a list of dge results, or list of lists with dge results
#'   in them. I sometimes iterate over each comparisons and have for each one
#'   a `dge` result, a `lfc` result (test against a threshold with treat),
#'   a `gsea` result, etc. If `dge.name` is provided, we assume `dge.res` is
#'   a list of lists and try to fish out the result with `dge.name`
#' @param fdr the FDR (padj) threshold for significance
#' @param lfc the minimum thrshold for significance. If you are summarizing
#'   a treat_lfc result, you will want to keep this at 0 because the threshold
#'   is already part of the test.
#' @param rank.down,rank.up Report the position of the genes listed here in the
#'   set of statistically significant down (or up) regulated genes. Use gene
#'   symbols for now.
#' @param dge.label the name of the colum in `comps` that stores the "label"
#'   of the result. We will assume `dge.res` is named by these values.
#' @examples
#' xs.all <- FacileData::an_fds() |>
#'   FacileData::samples() |> 
#'   FacileData::with_sample_covariates()
#' xs <- xs.all |> dplyr::filter(cond == "AKI")
#' comps <- dplyr::tibble(
#'   covariate = "cell_abbrev",
#'   numer = c("IMM", "EC"),
#'   denom = c("CNT", "IMM"),
#'   label = sprintf("%s_vs_%s", numer, denom)
#' )
#' fdge.all <- lapply(1:nrow(comps), function(i) {
#'   info <- comps[i,]
#'   xs |> 
#'     flm_def(info$covariate, numer = info$numer, denom = info$denom) |> 
#'     fdge(method = "voom", metadata = as.list(info))
#' })
#' names(fdge.all) <- comps$label
#'
#' fdsa <- summarize_fdge_all(comps, fdge.all, fdr = 0.1, lfc = log2(2))
#' rdsa <- summarize_dge_all(comps, rdge.all, fdr = 0.05, rank = c("HTT", "PMS1", "MSH3"))
summarize_fdge_all <- function(
  comps,
  dge_res,
  fdr = 0.1,
  lfc = 0,
  rank = NULL,
  ...,
  dge_name = "dge",
  dge_label = "label"
) {
  assert_data_frame(comps)
  assert_string(dge_label)
  labels <- assert_character(comps[[dge_label]])
  assert_list(dge_res, names = "unique")
  assert_subset(names(dge_res), labels)

  compstats <- dplyr::bind_cols(comps, {
    lapply(labels, function(l) {
      xres <- dge_res[[l]]
      if (!is(xres, "data.frame") && !is(xres, "FacileTtestAnalysisResult")) {
        xres <- xres[[dge_name]]
      }
      if (!is(xres, "data.frame") && !is(xres, "FacileTtestAnalysisResult")) {
        stop("Could not retrieve fdge result from embedded item named `", dge_name, "`")
        xres <- xres[[dge_name]]
      }
      summarize_fdge(xres, fdr = fdr, lfc = lfc, rank = rank, ...)
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_cols()
  compstats
}

#' Summarize a DGE result (data.frame) with dge counts and rank positions.
#'
#' For genes enumerated in `rank`, the user is given the option to
#' filter/include those based on their nominal pvalues, which is controlled
#' by the `rank_stats` option. By default, nominal pvalue is used for ranking.
#'
#' @export
#' @examples
#' xs.all <- FacileData::an_fds() |>
#'   FacileData::samples() |> 
#'   FacileData::with_sample_covariates()
#' xs <- xs.all |> dplyr::filter(cond == "AKI")
#' dge <- flm_def(xs, "cell_abbrev", "IMM", "EC") |> fdge(method = "voom")
#' summarize_fdge(dge, lfc = 1, fdr = 0.05)
#' (s1 <- summarize_fdge(dge, rank = c("LCP1", "SOX18")))
summarize_fdge <- function(
  x,
  fdr = 0.1,
  lfc = 1,
  rank = NULL,
  only_significant = TRUE,
  # rank_threshold = 0.05, rank_stats = c("pval", "padj"),
  ...,
  contrast = NULL
) {
  x. <- NULL
  if (is(x, "FacileTtestAnalysisResult")) {
    x. <- x
    x <- tidy(x)
  }
  assert_data_frame(x)
  assert_flag(only_significant)

  xs <- x |>
    dplyr::mutate(
      significant = (.data$padj <= .env$fdr) & (abs(.data$logFC) >= .env$lfc),
      keep = (.data$padj <= .env$fdr) & (abs(.data$logFC) >= .env$lfc)
    )

  rr <- NULL
  if (!is.null(rank)) {
    if (!is.data.frame(rank)) {
      assert_character(rank, min.len = 1)
      rr <- dplyr::filter(xs, .data$name %in% rank)
    } else {
      stopifnot("feature_id" %in% colnames(rank))
      rr <- dplyr::semi_join(xs, rank, by = "feature_id")
    }
  }

  if (!only_significant) {
    if (is.null(rank)) {
      stop("If only_significant = FALSE, then `rank` vector must be provided")
    }

    # down
    rrd <- dplyr::filter(rr, logFC < 0) |> dplyr::arrange(dplyr::desc(logFC))

    if (nrow(rrd) > 0L) {
      xs$keep <- xs$keep | (xs$logFC <= rrd$logFC[1] & xs$padj <= rrd$padj[1])
    }

    # up
    rru <- dplyr::filter(rr, logFC > 0) |> dplyr::arrange(logFC)
    if (nrow(rru) > 0L) {
      xs$keep <- xs$keep | (xs$logFC >= rru$logFC[1] & xs$padj <= rru$padj[1])
    }
  }

  xs <- xs |>
    dplyr::filter(.data$keep) |>
    dplyr::mutate(direction = ifelse(.data$logFC > 0, "up", "down"))
  if (nrow(xs) == 0L) {
    xs$direction <- character()
  }

  down <- xs |>
    dplyr::filter(direction == "down") |>
    dplyr::arrange(logFC)
  down <- dplyr::bind_rows(
    dplyr::filter(down, significant),
    dplyr::filter(down, !significant)
  ) |>
    dplyr::mutate(rank = -seq_len(n()))

  up <- xs |>
    dplyr::filter(direction == "up") |>
    dplyr::arrange(desc(logFC))
  up <- dplyr::bind_rows(
    dplyr::filter(up, significant),
    dplyr::filter(up, !significant)
  ) |>
    dplyr::mutate(rank = seq_len(n()))

  # stats <- tibble(
  #   dge_total = nrow(xs), dge_up = nrow(up), dge_down = nrow(down))
  stats <- tibble(
    dge_total = sum(xs$significant),
    dge_up = sum(xs$significant & xs$direction == "up"),
    dge_down = sum(xs$significant & xs$direction == "down")
  )

  ranked <- dplyr::bind_rows(up, down)

  if (!is.null(rr) && nrow(rr) > 0L) {
    xr <- rr |>
      dplyr::mutate(direction = ifelse(logFC > 0, "up", "down")) |>
      dplyr::left_join(select(ranked, feature_id, rank), by = "feature_id") |>
      dplyr::mutate(
        rank = ifelse(is.na(rank), sign(logFC) * Inf, rank)
      )

    # xr <- ranked |>
    #   dplyr::semi_join(rr, by = "feature_id") # dplyr::filter(ranked, .data$name %in% .env$rank | .data$feature_id %in% .env$rank)
    #
    # missed <- setdiff(rr$name, xr$name)
    # if (length(missed) > 0L) {
    #   missed <- dplyr::tibble(name = missed, logFC = 0, rank = NA_integer_)
    #   xr <- dplyr::bind_rows(xr, missed)
    # }

    sosmall <- 0.00001
    xr <- xr |>
      dplyr::mutate(
        # 99% means it's among the top 99% most regulated gene in a given direction
        rankpct = 1 -
          (abs(rank) / ifelse(rank < 0, stats$dge_down, stats$dge_up)),
        # rankpct is (-)Inf when there are 0 statistically sig DGE
        rankpct = ifelse(is.infinite(rankpct), sosmall, rankpct),
        rankpct = ifelse(
          # when we spanked on a gene because we were tracking/ranking it but it wasn't DGE
          # then rankpct will be < 0, and it will distort the rank.
          rankpct < 0 | is.na(rankpct) | rankpct == 0 | abs(rankpct) > 1,
          sosmall,
          rankpct
        ),
        rankpct = ifelse(logFC > 0, rankpct, -rankpct),
        rank = ifelse(is.na(rank), -Inf, rank),
        # pass_adj = .data$padj <= fdr & abs(logFC) >= lfc,
        pass_adj = .data$padj <= fdr,
        pass_adj = ifelse(is.na(pass_adj), FALSE, pass_adj),
        # pass_nom = .data$pval <= 0.05 & abs(logFC) >= lfc,
        pass_nom = .data$pval <= 0.05,
        pass_nom = ifelse(is.na(pass_nom), FALSE, pass_nom)
      )

    lfc.all <- c(
      xr$logFC,
      ranked$logFC[!ranked$name %in% xr$name],
      rnorm(nrow(x) - nrow(ranked), mean = 0, sd = 0.5)
    )
    # rnorm(nrow(x) - nrow(xs), mean = 0, sd = 0.5))
    # rep(0, 100))#nrow(x) - nrow(xs)))
    z.lfc <- scale(lfc.all, center = TRUE, scale = TRUE)[, 1]
    ranks <- xr |>
      dplyr::mutate(z = z.lfc[1:n()]) |>
      dplyr::select(name, z, rank, rankpct, logFC, pass_adj, pass_nom) |>
      tidyr::pivot_wider(
        names_from = "name",
        values_from = c(
          "z",
          "rank",
          "rankpct",
          "logFC",
          "pass_adj",
          "pass_nom"
        ),
        names_vary = "slowest"
      )
    names(ranks) <- paste0("dge_", names(ranks))

    # organize genes in the order they were passed in
    for (rname in rev(rank)) {
      ranks <- dplyr::relocate(ranks, dplyr::contains(rname), .before = 1L)
    }
    stats <- dplyr::bind_cols(stats, ranks)
  }

  stats
}

#' Stylized volcano plot(s)
#'
#' @export
#' @param x a `FacileTtestAnalysisResult` or a list of them to facet
#' @param lfc,fdrdefine significant result thresholds
#' @param highlight name or feature_id of genes to highlight. If the vector is
#'   named, then the names are the genes to highlight and the values are the
#'   colors, otherwise the values are the elements to highlight, and they will
#'   be highlighted using the color specified by `color.highlight`
#' @param label If `TRUE` (default), then we add labels to the highlighted
#'   genes using ggrepel
#' @param interactive If `TRUE` (default), the signficant and the highlighted
#'   genes will be available for "interactivity" using ggiraph. If the
#'   plot is plotted "normally", it is a standard ggplot, however you can call
#'   `ggiraph::girafe(ggobj = this)`, and you'll get an interactive plot
#' @param color.insig,size.insig,alpha.insig aesthetics of points of insig genes
#' @param color.sig,size.sig,alpha.sig aesthetics of points of significant genes
#' @param color.highlight,size.highlight aesthetics of genes to highlight
#' @param xlims if `"square"` (default), will set x axis to be max(abs(logFC))
#' @param ... stuff
#' @examples
#' xs.all <- FacileData::an_fds() |>
#'   FacileData::samples() |> 
#'   FacileData::with_sample_covariates()
#' xs <- xs.all |> dplyr::filter(cond == "AKI")
#' dge <- flm_def(xs, "cell_abbrev", "IMM", "EC") |> fdge(method = "voom")
#' 
#' res <- tidy(dge) |> dplyr::arrange(pval)
#' hi <- res |> arrange(desc(logFC)) |> head(100) |> sample_n(20)
#' lbl <- head(res, 5)
#' 
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = "ENSG00000067113"
#' )
#' 
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.05,
#'   highlight = hi,
#'   label = lbl
#' )
#'
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = c(PLPP1 = "red", ROBO4 = "pink", LCP1 = "orange")
#' )
#' 
#' 
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = head(res, 5)
#' )
#' 
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = head(res, 5),
#'   label = FALSE
#' )
#'
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = head(res, 5),
#'   label = TRUE
#' )
#' 
#' fvolcano(
#'   dge,
#'   lfc = log2(2.5),
#'   fdr = 0.1,
#'   highlight = head(res, 5),
#'   label = "HCST"
#' )
fvolcano <- function(
  x,
  lfc = 1,
  fdr = 0.10,
  highlight = NULL,
  label = !is.null(highlight),
  label_max_overlaps = 20,
  with_stats = FALSE,
  nudge_x = 0,
  nudge_y = 0,
  interactive = FALSE,
  interactive_opts = list(),
  color.insig = "lightgrey",
  size.insig = 1,
  alpha.insig = 0.7,
  color.sig = "cornflowerblue",
  size.sig = 1,
  alpha.sig = 0.7,
  color.highlight = "firebrick2",
  size.highlight = 2.5,
  color.label = color.highlight,
  size.label = size.highlight,
  size.label_text = 4,
  xlims = "square",
  ylims = NULL,
  vline = NULL,
  vline_color = "firebrick",
  vline_type = "dashed",
  event_source = "A",
  ...,
  facet_ncol = 5,
  title = NULL,
  yaxis = c("padj", "pval"),
  contrast = NULL,
  dge_method = "DESeq2",
  grouping_var = "portfolio_id"
) {
  yaxis <- match.arg(yaxis)

  # do we want to separate the label from the highlight?
  callout.label <- color.highlight != color.label

  checkmate::assert_flag(interactive)
  if (is(x, "data.frame") || is(x, "FacileTtestAnalysisResult")) {
    x <- list(dge = x)
  }
  x <- lapply(x, function(y) {
    if (is(y, "FacileTtestAnalysisResult")) tidy(y) else y
  })
  stopifnot(all(sapply(x, test_data_frame)))

  vdat.all <- lapply(names(x), function(xname) {
    out <- x[[xname]]
    out[["rname"]] <- xname
    out[["grouping"]] <- out[[grouping_var]][1L]
    out
  }) |>
    dplyr::bind_rows()

  # Depending on where we pull the stats from, the gene names are stored in
  # "symbol", "name", or "gene" columns. Let's see what is the what
  gname <- "name"
  if (!"name" %in% colnames(vdat.all)) {
    gname <- if ("symbol" %in% colnames(vdat.all)) "symbol" else "gene"
  }
  if (!is.character(vdat.all[[gname]])) {
    stop("Can't find 'gene name' column in these stats tables")
  }
  vdat.all$feature_name <- vdat.all[[gname]]
  if (!is.numeric(vdat.all$pval)) {
    vdat.all$pval <- 1
  }
  if (!is.numeric(vdat.all$padj)) {
    vdat.all$padj <- 1
  }

  vdat.all <- vdat.all |>
    dplyr::mutate(
      rname = factor(rname, unique(rname)),
      sig = .data$padj <= .env$fdr & abs(.data$logFC) >= lfc
    ) |>
    dplyr::select(
      rname,
      # {{ gname }},
      feature_id,
      feature_name,
      logFC,
      pval,
      padj,
      sig,
      # highlight,
      everything()
    ) |>
    dplyr::mutate(
      .group = ifelse(sig, "sig", "insig"),
      tooltip = sprintf(
        "%s\nlogFC: %0.2f\npadj: %0.2f\npval: %0.2f",
        feature_name,
        logFC,
        padj,
        pval
      ),
      label = ""
    )

  if (is.null(highlight) && !checkmate::test_flag(label)) {
    highlight <- label
    label <- TRUE
  }
  
  if (is.data.frame(highlight)) {
    if (is.character(highlight[["feature_id"]])) {
      h <- highlight[["feature_id"]]
    } else {
      h <- highlight[[gname]]
    }
    if (is.null(h)) {
      warning("Can't figure out what to do with `highlight`, ignoring.")
      highlight <- NULL
    } else {
      if (is.null(highlight[["color"]])) {
        highlight <- h
      } else {
        highlight <- setNames(highlight[["color"]], h)
      }
    }
  }
  
  # Change `highlight` vector to have names() of targets, and values as color
  # User can pass either gene name/symbol or feature_id
  hfeatures <- NULL
  if (test_character(highlight, min.len = 1L)) {
    if (is.null(names(highlight))) {
      # user just passed in names of things to highlight
      highlight <- setNames(rep(color.highlight, length(highlight)), highlight)
    }
    use.def.color <- nchar(names(highlight)) == 0L
    if (any(use.def.color)) {
      highlight[use.def.color] <- color.highlight
    }
    highlight <- highlight[!duplicated(names(highlight))]
    assert_character(highlight, names = "unique")
    hfeatures <- feature_lookup(
      names(highlight),
      vdat.all,
      feature_name = "feature_name"
    )
    if (nrow(hfeatures)) {
      vdat.all$.group[hfeatures$feature_row_index] <- hfeatures$feature_query
    }
  }
  
  # anything to label?
  if (isTRUE(label)) {
    label <- hfeatures
  } else if (!isFALSE(label)) {
    lbl <- NULL
    if (is.data.frame(label)) {
      lbl <- label[[gname]]
      if (is.null(lbl)) {
        lbl <- label[["feature_id"]]
      }
    } else if (is.character(label)) {
      lbl <- label
    }
    if (is.null(lbl)) {
      warning("Don't know how to handle `label`, skipping")
      label <- FALSE
    } else {
      label <- feature_lookup(
        lbl,
        vdat.all,
        feature_name = "feature_name"
      )
    }
  }
  if (is.data.frame(label)) {
    if (with_stats) {
      vdat.all$label[label$feature_row_index] <- vdat.all$tooltip[label$feature_row_index]  
    } else {
      vdat.all$label[label$feature_row_index] <- vdat.all$feature_name[label$feature_row_index]  
    }
    label <- TRUE
  }
  
  vdat <- dplyr::bind_rows(
    dplyr::filter(vdat.all, .group == "insig"),
    dplyr::filter(vdat.all, .group == "sig"),
    dplyr::filter(vdat.all, !.group %in% c("insig", "sig"))
  ) |>
    dplyr::mutate(
      data_id = seq_along(label)
    )
  
  if (yaxis == "padj") {
    vdat$yax <- -log10(vdat$padj + 1e-6)
  } else {
    vdat$yax <- -log10(vdat$pval + 1e-6)
  }

  ggcolors <- c(insig = color.insig, sig = color.sig)
  ggcolors <- c(ggcolors, highlight)
  ggcolors <- ggcolors[!duplicated(names(ggcolors))]

  ggalpha <- c(insig = alpha.insig, sig = alpha.sig)
  ggalpha <- c(ggalpha, setNames(rep(1, length(highlight)), names(highlight)))
  ggalpha <- ggalpha[!duplicated(names(ggalpha))]

  ggsize <- c(insig = size.insig, sig = size.sig)
  ggsize <- c(
    ggsize,
    setNames(rep(size.highlight, length(highlight)), names(highlight))
  )
  ggsize <- ggsize[!duplicated(names(ggsize))]

  gg <- ggplot2::ggplot(vdat) +
    ggplot2::aes(x = logFC, y = yax) +
    ggplot2::labs(
      y = sprintf("-log10(%s)", if (yaxis == "pval") "p-value" else "FDR"),
      x = "log2 fold change"
    )

  if (length(x) > 1L) {
    gg <- gg + ggplot2::facet_wrap(~rname, ncol = facet_ncol)
  }

  if (interactive) {
    gg <- gg +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .group,
          alpha = .group,
          size = .group
        ),
        data = dplyr::filter(vdat, .group == "insig")
      ) +
      ggiraph::geom_point_interactive(
        ggplot2::aes(
          color = .group,
          alpha = .group,
          size = .group,
          tooltip = tooltip,
          data_id = feature_id
        ),
        data = dplyr::filter(vdat, .group != "insig")
      )
  } else {
    gg <- gg +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .group,
          alpha = .group,
          size = .group
        )
      )
  }

  if (checkmate::test_number(vline)) {
    gg <- gg +
      geom_vline(
        xintercept = vline,
        color = vline_color,
        linetype = vline_type
      )
  }
  if (label) {
    lbl.dat <- dplyr::filter(vdat, label != "")
    gg <- gg +
      ggplot2::geom_point(
        data = lbl.dat,
        color = color.label,
        size = size.label
      ) +
      ggrepel::geom_label_repel(
        ggplot2::aes(label = label),
        max.overlaps = label_max_overlaps,
        nudge_x = nudge_x,
        nudge_y = nudge_y,
        size = size.label_text,
        # color = color.label,
        data = lbl.dat
      )
  }

  gg <- gg +
    ggplot2::scale_color_manual(values = ggcolors) +
    ggplot2::scale_alpha_manual(values = ggalpha) +
    ggplot2::scale_size_manual(values = ggsize)

  if (test_string(xlims) && xlims == "square") {
    # xlims <- max(ceiling(abs(vdat$logFC) + 0.2))
    xlims <- max(abs(vdat$logFC)) + 0.1
    xlims <- c(-xlims, xlims)
  }
  if (is.numeric(xlims)) {
    if (length(xlims) != 2) {
      warning("invalid value for xlims -- must be numeric vector of length 2")
    } else {
      breaks <- seq(floor(xlims[1]), ceiling(xlims[2]), by = 1)
      gg <- gg +
        ggplot2::scale_x_continuous(
          limits = xlims,
          breaks = breaks
        )
    }
  }
  if (!test_numeric(ylims, len = 2)) {
    ylims <- c(0, max(vdat$yax))
  }
  gg <- gg + ggplot2::scale_y_continuous(limits = ylims)

  gg <- gg + ggplot2::theme(legend.position = "none")
  if (test_string(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }
  attr(gg, "data") <- vdat
  gg
}
