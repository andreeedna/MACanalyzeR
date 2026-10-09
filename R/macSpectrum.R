## MacSpectrum / AtheroSpectrum indices (MPI, AMDI, MDFI) for every cell.
.macspec <- function(mac_obj) {
  org <- mac_obj@organism
  refs <- list(M1 = MACanalyzeR::M1_mean_new, M2 = MACanalyzeR::M2_mean_new, M0 = MACanalyzeR::M0_mean_new,
               con = MACanalyzeR::con_mean_new, foam = MACanalyzeR::foam_mean_new)

  ## reference profiles are defined on mouse genes
  if (org != "mm") {
    dict <- MACanalyzeR::human_mouse_dict
    refs <- lapply(refs, function(r) {
      names(r) <- dict[names(r), "hs_symbol"]
      r[!is.na(names(r))]
    })
  }

  genes <- unique(unlist(lapply(refs, names)))
  mac_mtx <- .expr(mac_obj, genes)
  mac_mtx <- mac_mtx - rowMeans(mac_mtx)

  ## correlation-like similarity between each cell and a reference profile
  ## (genes centered across cells, root-mean-square scaling as in MacSpectrum)
  sim <- function(ref_name, gene_ref = ref_name) {
    g <- intersect(intersect(names(refs[[gene_ref]]), names(refs[[ref_name]])), rownames(mac_mtx))
    if (length(g) < 10) stop("Too few reference genes found in the data for MacSpectrum (", length(g), ")", call. = FALSE)
    mtx <- mac_mtx[g, , drop = FALSE]
    ref <- refs[[ref_name]][g]
    n <- length(g)
    cell_sigma <- sqrt(colSums(mtx^2) / n)
    ref_sigma <- sqrt(sum(ref^2) / n)
    colSums(sweep(mtx, 2, cell_sigma, "/") * (ref / ref_sigma)) / n
  }

  ## projection of (x0, y0) on the line a*x + b*y + c = 0, rescaled to [-50, 50]
  project <- function(x0, y0, a, b, c) {
    d_sqr <- (a * x0 + b * y0 + c)^2 / (a^2 + b^2)
    x_start <- -1; y_start <- (-a) * x_start + (-c)
    x_end <- 1;    y_end <- (-a) * x_end + (-c)
    l <- ((x0 - x_start)^2 + (y0 - y_start)^2 - d_sqr)^0.5
    l_max <- ((x_end - x_start)^2 + (y_end - y_start)^2 - d_sqr)^0.5
    l / l_max * 100 - 50
  }

  M1_pearson <- sim("M1")
  M2_pearson <- sim("M2", gene_ref = "M1")
  M0_pearson <- sim("M0")
  foam_pearson <- sim("foam")
  con_pearson <- sim("con", gene_ref = "foam")

  data.frame(MPI = project(M1_pearson, M2_pearson, a = 0.991414467, b = 1, c = -0.0185412856),
             AMDI = -M0_pearson * 50,
             MDFI = project(foam_pearson, con_pearson, a = 1, b = 1, c = 0),
             row.names = colnames(mac_mtx))
}


MacSpectrumIndex <- function(mac_obj) {
  .check_obj(mac_obj)
  ms <- .macspec(mac_obj)
  for (v in colnames(ms)) mac_obj@MetaData[[v]] <- ms[[v]]
  mac_obj
}


macSpectrum <- function(mac_obj, mode = "macspectrum", plot.by = mac_obj@ident, split.by = NULL, style = "density",
                        txt.size = 17, col = NULL, ncol = 3) {
  .check_obj(mac_obj)
  if (!style %in% c("point", "density")) {
    stop("The style argument must be 'point' (scatter plot) or 'density' (2D density plot)", call. = FALSE)
  }
  if (!mode %in% c("macspectrum", "atherospectrum")) {
    stop("The mode argument must be 'macspectrum' or 'atherospectrum'", call. = FALSE)
  }
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")

  dfplot <- mac_obj@MetaData
  if (!all(c("MPI", "AMDI", "MDFI") %in% colnames(dfplot))) {
    ms <- .macspec(mac_obj)
    for (v in colnames(ms)) dfplot[[v]] <- ms[[v]]
  }

  yvar <- if (mode == "macspectrum") "AMDI" else "MDFI"
  plot <- ggplot(dfplot, aes(x = .data[["MPI"]], y = .data[[yvar]])) +
    labs(fill = plot.by, color = plot.by) +
    geom_vline(xintercept = 1, colour = "grey30") +
    geom_hline(yintercept = 0.5, colour = "grey30") +
    .mac_theme(txt.size)

  if (style == "point") {
    plot <- plot + geom_point(aes(color = .data[[plot.by]]), size = 2)
  } else {
    plot <- plot + stat_density_2d(aes(color = .data[[plot.by]], fill = .data[[plot.by]]),
                                   geom = "polygon", alpha = 0.2, position = "identity")
  }

  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)

  if (!is.null(col)) {
    if (length(col) == length(unique(stats::na.omit(dfplot[[plot.by]])))) {
      plot <- plot + scale_color_manual(values = col) + scale_fill_manual(values = col)
    } else {
      warning("The number of colours must match the number of groups of ", plot.by, ": using the default palette",
              call. = FALSE)
    }
  }
  plot
}
