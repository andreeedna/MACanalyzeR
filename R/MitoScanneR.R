#' Mitochondrial Scanning for Macrophages
#'
#' @param object A Seurat object
#' @param group.by Metadata column for grouping (default: from CreateMacObj)
#' @param assay Assay to use (default: DefaultAssay(object))
#' @param layer Layer to use (default: "data")
#' @return A Seurat object with mitochondrial scores
#' @importFrom stats weighted.mean var
#' @importFrom future.apply future_lapply
#' @importFrom Seurat GetAssayData DefaultAssay
#' @export
MitoScanneR <- function(object, group.by = NULL, assay = NULL, layer = "data") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  if (is.null(group.by)) {
    group.by <- mac_info$ident
  }
  if (is.null(object@meta.data[[group.by]])) {
    stop("Metadata column '", group.by, "' not found")
  }

  if (mac_info$organism == "mm") {
    MitoPath <- MACanalyzeR::mouse_MitoCarta
  } else {
    MitoPath <- MACanalyzeR::human_MitoCarta
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  cond_tot <- object@meta.data[[group.by]]
  cond <- unique(cond_tot)

  check.cols <- function(vect) {
    all(vect >= 0.001)
  }

  score_mito <- function(p_name, category_pathways) {
    genes <- category_pathways[[p_name]]
    genes <- intersect(genes, rownames(mtx))
    if (length(genes) < 5) {
      return(NULL)
    }

    mtx_pathway <- mtx[genes, , drop = FALSE]
    mtx_pathway <- mtx_pathway[Matrix::rowSums(mtx_pathway) > 0, , drop = FALSE]

    mtx_dense <- as.matrix(mtx_pathway)
    sample_mean <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if (length(keep) < 3) {
      return(NULL)
    }

    mtx_dense <- mtx_dense[keep, , drop = FALSE]
    mean_expr_by_cluster <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    ratio_expr <- t(mean_expr_by_cluster) / colMeans(mean_expr_by_cluster)
    gene_weight <- apply(mtx_dense, 1, stats::var)
    mean_exp_pathway <- apply(ratio_expr, 2, function(x) stats::weighted.mean(x, gene_weight / sum(gene_weight)))
    return(mean_exp_pathway)
  }

  MitoList <- list()
  for (m in names(MitoPath)) {
    message("Scanning category: ", m)
    category_pathways <- MitoPath[[m]]
    path_names <- names(category_pathways)

    scored_list <- future.apply::future_lapply(path_names, score_mito, category_pathways = category_pathways, future.seed = TRUE)
    names(scored_list) <- path_names
    scored_list <- scored_list[!sapply(scored_list, is.null)]

    MitoList[[m]] <- as.data.frame(do.call(rbind, scored_list))
  }

  object@misc$MACanalyzeR$MitoScanneR[[group.by]] <- MitoList
  return(object)
}

#' Heatmap of Mitochondrial Scores
#'
#' @param object A Seurat object
#' @param group.by Metadata column used (default: from CreateMacObj)
#' @param cols Colors for heatmap
#' @param base.size Font size
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
HeatMito <- function(object, group.by = NULL, cols = c("blue", "white", "red"), base.size = 10) {
  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  mito_data <- object@misc$MACanalyzeR$MitoScanneR[[group.by]]
  if (is.null(mito_data)) {
    stop("Run 'MitoScanneR' first")
  }

  for (mitoname in names(mito_data)[length(mito_data):1]) {
    dat <- mito_data[[mitoname]]
    sort_row <- c()
    for (i in colnames(dat)) {
      select_row <- which(apply(dat, 1, max) == dat[, i])
      tmp <- rownames(dat)[select_row][order(dat[select_row, i], decreasing = TRUE)]
      sort_row <- c(sort_row, tmp)
    }
    dat[is.na(dat)] <- 1
    b <- max((max(dat) - 1), (1 - min(dat)))
    mybreaks <- c(
      seq((1 - b), 1 - (b / 3), length.out = 33),
      seq(1 - (b / 3) + 0.01, 1 + (b / 3), length.out = 34),
      seq(1 + (b / 3) + 0.01, 1 + b, length.out = 33)
    )
    color <- grDevices::colorRampPalette(cols)(100)
    pheatmap::pheatmap(dat[sort_row, ], cluster_cols = FALSE, cluster_rows = FALSE, color = color, breaks = mybreaks, main = gsub("_", " ", mitoname),
             fontsize = base.size, border_color = FALSE)
  }
}

#' Mitochondrial/Nuclear Balance
#'
#' @param object A Seurat object
#' @param group.by Metadata column for grouping
#' @param base.size Font size
#' @param ncol Number of columns for facetting
#' @param cols Custom colors
#' @param intercept Show intercept at 0
#' @param title Plot title
#' @param assay Assay (default: DefaultAssay(object))
#' @param layer Layer
#' @return A ggplot object
#' @import ggplot2
#' @importFrom Seurat DefaultAssay
#' @export
MitoBalance <- function(object, group.by = "Sample", base.size = 15, ncol = 1, cols = NULL, intercept = FALSE, title = "MitoNuclear Balance", assay = NULL, layer = "data") {
  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  cond_tot <- as.vector(object@meta.data[[group.by]])
  path_names <- c("OXPHOS mitochondrial subunits", "OXPHOS nuclear subunits")

  if (mac_info$organism == "mm") {
    MitoPath <- MACanalyzeR::mouse_MitoCarta$OXPHOS
  } else {
    MitoPath <- MACanalyzeR::human_MitoCarta$OXPHOS
  }

  check.cols <- function(vect) {
    all(vect >= 0.001)
  }

  calc_score <- function(p) {
    genes <- intersect(MitoPath[[p]], rownames(mtx))
    if (length(genes) < 5) {
      return(NULL)
    }
    mtx_p <- mtx[genes, , drop = FALSE]
    mtx_p <- mtx_p[Matrix::rowSums(mtx_p) > 0, , drop = FALSE]

    mtx_dense <- as.matrix(mtx_p)
    sample_mean <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if (length(keep) < 3) {
      return(NULL)
    }

    mtx_dense <- mtx_dense[keep, , drop = FALSE]
    mean_expr <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    ratio_expr <- mtx_dense / colMeans(mean_expr)
    gene_weight <- apply(mtx_dense, 1, stats::var)
    return(apply(ratio_expr, 2, function(x) stats::weighted.mean(x, gene_weight / sum(gene_weight))))
  }

  scores <- lapply(path_names, calc_score)
  if (any(sapply(scores, is.null))) {
    stop("Could not calculate scores for both subunits")
  }

  dat <- data.frame(Mito = scores[[1]], Nuclear = scores[[2]])
  dat$MitochondrialBalance <- log10(dat$Mito / dat$Nuclear)
  dat[[group.by]] <- object@meta.data[[group.by]]

  threshold <- is.finite(dat$MitochondrialBalance)
  if (any(!threshold)) {
    message("Removed ", sum(!threshold), " cells containing non-finite values")
  }
  dat <- dat[threshold, ]

  p1 <- ggplot(dat, aes(x = .data[[group.by]], y = MitochondrialBalance, fill = .data[[group.by]], color = .data[[group.by]])) +
    geom_violin() +
    geom_boxplot(width = 0.3, color = "grey30", alpha = 0.2) +
    labs(x = "", y = "log10(Mitochondrial/Nuclear GeneRatio)") +
    ggtitle(title) +
    theme_classic(base_size = base.size) +
    theme(legend.position = "none")

  if (intercept) {
    p1 <- p1 + geom_hline(yintercept = 0, colour = 'grey30')
  }
  if (!is.null(cols)) {
    p1 <- p1 + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
  }

  return(p1)
}

#' Glycolysis/OXPHOS Balance
#'
#' @param object A Seurat object
#' @param group.by Metadata column for grouping
#' @param base.size Font size
#' @param ncol Number of columns for facetting
#' @param cols Custom colors
#' @param intercept Show intercept at 0
#' @param title Plot title
#' @param assay Assay (default: DefaultAssay(object))
#' @param layer Layer
#' @return A ggplot object
#' @import ggplot2
#' @importFrom Seurat DefaultAssay
#' @export
GliOxBalance <- function(object, group.by = "Cluster", base.size = 15, ncol = 1, cols = NULL, intercept = FALSE, title = "GlycoOxphos Balance", assay = NULL, layer = "data") {
  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  cond_tot <- as.vector(object@meta.data[[group.by]])
  path_names <- c("GLYCOLYSIS_GLUCONEOGENESIS", "OXIDATIVE_PHOSPHORYLATION")

  if (mac_info$organism == "mm") {
    PathwaySet <- MACanalyzeR::mouse_selected
  } else {
    PathwaySet <- MACanalyzeR::human_selected
  }

  check.cols <- function(vect) {
    all(vect >= 0.001)
  }

  calc_score <- function(p) {
    genes <- intersect(PathwaySet[[p]], rownames(mtx))
    if (length(genes) < 5) {
      return(NULL)
    }
    mtx_p <- mtx[genes, , drop = FALSE]
    mtx_p <- mtx_p[Matrix::rowSums(mtx_p) > 0, , drop = FALSE]

    mtx_dense <- as.matrix(mtx_p)
    sample_mean <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if (length(keep) < 3) {
      return(NULL)
    }

    mtx_dense <- mtx_dense[keep, , drop = FALSE]
    mean_expr <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))
    ratio_expr <- mtx_dense / colMeans(mean_expr)
    gene_weight <- apply(mtx_dense, 1, stats::var)
    return(apply(ratio_expr, 2, function(x) stats::weighted.mean(x, gene_weight / sum(gene_weight))))
  }

  scores <- lapply(path_names, calc_score)
  if (any(sapply(scores, is.null))) {
    stop("Could not calculate scores for both pathways")
  }

  dat <- data.frame(Glyco = scores[[1]], OXPHOS = scores[[2]])
  dat$Balance <- log10(dat$Glyco / dat$OXPHOS)
  dat[[group.by]] <- object@meta.data[[group.by]]

  threshold <- is.finite(dat$Balance)
  if (any(!threshold)) {
    message("Removed ", sum(!threshold), " cells containing non-finite values")
  }
  dat <- dat[threshold, ]

  p1 <- ggplot(dat, aes(x = .data[[group.by]], y = Balance, fill = .data[[group.by]], color = .data[[group.by]])) +
    geom_violin() +
    geom_boxplot(width = 0.3, color = "grey30", alpha = 0.2) +
    labs(x = "", y = "log10(Glycolysis/OXPHOS GeneRatio)") +
    ggtitle(title) +
    theme_classic(base_size = base.size) +
    theme(legend.position = "none")

  if (intercept) {
    p1 <- p1 + geom_hline(yintercept = 0, colour = 'grey30')
  }
  if (!is.null(cols)) {
    p1 <- p1 + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
  }

  return(p1)
}

