MitoScanneR <- function(mac_obj, meta = mac_obj@ident, verbose = TRUE) {
  .check_obj(mac_obj)
  .check_col(mac_obj, meta, "meta")

  MitoPath <- .org_data(mac_obj, MACanalyzeR::mouse_MitoCarta, MACanalyzeR::human_MitoCarta)
  groups <- mac_obj@MetaData[[meta]]

  MitoList <- list()
  for (m in names(MitoPath)) {
    res <- .score_list(mac_obj, MitoPath[[m]], groups, single = FALSE, test = FALSE,
                       verbose = verbose, label = m)
    MitoList[[m]] <- as.data.frame(res$total)
  }
  mac_obj@MitoScanneR[[meta]] <- MitoList
  mac_obj
}

HeatMito <- function(mac_obj, plot.by = mac_obj@ident, col = c("blue", "white", "red"), txt.size = 10) {
  .check_obj(mac_obj)
  res <- mac_obj@MitoScanneR[[plot.by]]
  if (length(res) == 0) {
    stop("No MitoScanneR results for '", plot.by, "'. Run MitoScanneR(mac_obj, meta = '", plot.by, "') first.",
         call. = FALSE)
  }

  for (mitoname in rev(names(res))) {
    dat <- as.matrix(res[[mitoname]])
    dat <- dat[rowSums(!is.na(dat)) > 0, , drop = FALSE]
    if (nrow(dat) == 0) {
      message("No pathway of '", mitoname, "' could be scored: skipped")
      next
    }
    dat <- dat[.heat_order(dat), , drop = FALSE]
    dat[is.na(dat)] <- 1
    pheatmap(dat, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(col)(100),
             breaks = .heat_breaks(dat), main = gsub("_", " ", mitoname), fontsize = txt.size, border_color = NA)
  }
  invisible(NULL)
}


## Per-cell balance between two gene sets: log10(score set1 / score set2).
.balance_data <- function(mac_obj, set1, set2, plot.by, names = c("set1", "set2")) {
  groups <- mac_obj@MetaData[[plot.by]]
  sc <- lapply(list(set1, set2), function(g) {
    s <- .score_set(.expr(mac_obj, g), groups, single = TRUE)
    if (is.null(s)) NULL else s$single
  })
  for (i in 1:2) {
    if (is.null(sc[[i]])) {
      found <- length(intersect(list(set1, set2)[[i]], rownames(mac_obj@CountData@NormalizedCount)))
      stop("The gene set '", names[i], "' could not be scored: ", found, " of ", length(list(set1, set2)[[i]]),
           " genes found in the data (at least 5 detected genes, and 3 expressed in every group, are needed)",
           call. = FALSE)
    }
  }
  dat <- mac_obj@MetaData
  dat$.set1 <- sc[[1]]
  dat$.set2 <- sc[[2]]
  dat$.balance <- log10(sc[[1]] / sc[[2]])
  ok <- is.finite(dat$.balance)
  if (sum(!ok) > 0) message("Removed ", sum(!ok), " cells with a non-finite value (score equal to 0)")
  dat[ok, , drop = FALSE]
}


GeneBalance <- function(mac_obj, set1, set2, plot.by = mac_obj@ident, split.by = NULL, txt.size = 15,
                        ncol = 3, col = NULL, intercept = FALSE, title = NULL,
                        ylab = "log10(set1 / set2) score ratio") {
  .check_obj(mac_obj)
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")
  if (!is.character(set1) || !is.character(set2)) stop("'set1' and 'set2' must be vectors of gene names", call. = FALSE)

  dat <- .balance_data(mac_obj, set1, set2, plot.by, names = c("set1", "set2"))
  dat[[plot.by]] <- factor(dat[[plot.by]])
  dat <- dat[!is.na(dat[[plot.by]]), , drop = FALSE]

  p1 <- ggplot(dat, aes(x = .data[[plot.by]], y = .data[[".balance"]],
                        fill = .data[[plot.by]], color = .data[[plot.by]])) +
    geom_violin() +
    geom_boxplot(width = 0.3, color = "grey30", alpha = 0.2) +
    labs(x = "", y = ylab) +
    .mac_theme(txt.size) +
    theme(legend.position = "none")

  if (!is.null(title)) p1 <- p1 + ggtitle(title)
  if (isTRUE(intercept)) p1 <- p1 + geom_hline(yintercept = 0, colour = "grey30")
  if (.check_colors(col, dat[[plot.by]], "plot.by")) {
    p1 <- p1 + scale_fill_manual(values = col) + scale_color_manual(values = col)
  }
  if (!is.null(split.by)) p1 <- p1 + facet_wrap(split.by, ncol = ncol)
  p1
}


MitoBalance <- function(mac_obj, plot.by = mac_obj@ident, txt.size = 15, ncol = 3, col = NULL,
                        intercept = FALSE, title = "Mito-nuclear balance", split.by = NULL) {
  .check_obj(mac_obj)
  oxphos <- .org_data(mac_obj, MACanalyzeR::mouse_MitoCarta, MACanalyzeR::human_MitoCarta)$OXPHOS
  mito <- oxphos[["OXPHOS mitochondrial subunits"]]
  nuc <- oxphos[["OXPHOS nuclear subunits"]]

  n_mito <- length(intersect(mito, rownames(mac_obj@CountData@NormalizedCount)))
  if (n_mito < 5) {
    stop("Only ", n_mito, " of the ", length(mito), " mitochondrially encoded OXPHOS genes (",
         if (identical(mac_obj@organism, "hs")) "MT-" else "mt-",
         " genes) are present in the data. They are often removed during quality control: ",
         "keep them in the Seurat object to use MitoBalance().", call. = FALSE)
  }

  GeneBalance(mac_obj, mito, nuc, plot.by = plot.by, split.by = split.by, txt.size = txt.size, ncol = ncol,
              col = col, intercept = intercept, title = title,
              ylab = "log10(mitochondrial / nuclear OXPHOS)")
}


GliOxBalance <- function(mac_obj, plot.by = mac_obj@ident, txt.size = 15, ncol = 3, col = NULL,
                         intercept = FALSE, title = "Glycolysis/OXPHOS balance", split.by = NULL) {
  .check_obj(mac_obj)
  sel <- .org_data(mac_obj, MACanalyzeR::mouse_selected, MACanalyzeR::human_selected)

  GeneBalance(mac_obj, sel[["GLYCOLYSIS_GLUCONEOGENESIS"]], sel[["OXIDATIVE_PHOSPHORYLATION"]],
              plot.by = plot.by, split.by = split.by, txt.size = txt.size, ncol = ncol, col = col,
              intercept = intercept, title = title, ylab = "log10(glycolysis / OXPHOS)")
}
