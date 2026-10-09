#' Macrophage Polarization Analysis
#'
#' @param object A Seurat object
#' @param assay Assay to use (default: DefaultAssay(object))
#' @param layer Layer to use (default: "data")
#' @return A Seurat object with polarization metadata
#' @importFrom stats kmeans weighted.mean var
#' @importFrom Seurat GetAssayData DefaultAssay
#' @export
MacPolarizeR <- function(object, assay = NULL, layer = "data") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  organism <- mac_info$organism
  assay <- assay %||% Seurat::DefaultAssay(object)

  if (organism == "mm") {
    inf <- MACanalyzeR::mouse_PolGenes
  } else {
    inf <- MACanalyzeR::human_PolGenes
  }

  mac_mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  m1 <- intersect(inf$M1, rownames(mac_mtx))
  m2 <- intersect(inf$M2, rownames(mac_mtx))
  genes <- c(m1, m2)

  # Subset and convert to dense for kmeans (only for the few polarization genes)
  m_mtx <- t(as.matrix(mac_mtx[genes, ]))
  mac_km <- stats::kmeans(m_mtx, 3, iter.max = 5000)

  polarization_results <- list()

  for (store in c("Total", "SingleCell")) {
    if (store == "SingleCell") {
      path_expression <- matrix(NA, nrow = ncol(mac_mtx), ncol = 2,
                                dimnames = list(colnames(mac_mtx), c("M1", "M2")))
    } else {
      path_expression <- matrix(NA, nrow = 2, ncol = 3,
                                dimnames = list(c("M1", "M2"), c("1", "2", "3")))
    }

    for (p in c("M1", "M2")) {
      p_genes <- intersect(inf[[p]], rownames(mac_mtx))
      mtx_pathway <- mac_mtx[p_genes, , drop = FALSE]

      # Use rowSums on sparse matrix
      mtx_pathway <- mtx_pathway[Matrix::rowSums(mtx_pathway) > 0, , drop = FALSE]

      # Optimization: calculate cluster means efficiently
      # Using a dense subset for these calculations as p_genes is usually small
      mtx_dense <- as.matrix(mtx_pathway)

      if (store == "SingleCell") {
        mean_expr_by_cluster <- apply(mtx_dense, 1, function(x) by(x, mac_km$cluster, mean))
        # mean_expr_by_cluster is clusters x genes
        # colMeans(mean_expr_by_cluster) is mean of cluster-means for each gene
        ratio_expr <- mtx_dense / colMeans(mean_expr_by_cluster)
        gene_weight <- apply(mtx_dense, 1, stats::var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) stats::weighted.mean(x, gene_weight / sum(gene_weight)))
        path_expression[, p] <- mean_exp_pathway
      } else {
        sample_mean <- apply(mtx_dense, 1, function(x) by(x, mac_km$cluster, mean))
        ratio_expr <- t(sample_mean) / colMeans(sample_mean)
        gene_weight <- apply(mtx_dense, 1, stats::var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) stats::weighted.mean(x, gene_weight / sum(gene_weight)))
        path_expression[p, ] <- mean_exp_pathway
      }
    }
    polarization_results[[store]] <- path_expression
  }

  mi <- which.max(polarization_results$Total["M1", ])
  li <- which.max(polarization_results$Total["M2", ])

  cluster_names <- c("1", "2", "3")
  if (mi == li) {
    cluster_names[mi] <- "Inflammatory"
    remaining <- setdiff(1:3, mi)
    li_new <- remaining[which.max(polarization_results$Total["M2", remaining])]
    cluster_names[li_new] <- "Healing"
    cluster_names[setdiff(1:3, c(mi, li_new))] <- "Transitional"
  } else {
    cluster_names[mi] <- "Inflammatory"
    cluster_names[li] <- "Healing"
    cluster_names[setdiff(1:3, c(mi, li))] <- "Transitional"
  }

  colnames(polarization_results$Total) <- cluster_names

  object$MacPolarizeR <- factor(mac_km$cluster, levels = 1:3, labels = cluster_names)
  object$MacPolarizeR <- factor(object$MacPolarizeR, levels = c("Inflammatory", "Transitional", "Healing"))

  object@misc$MACanalyzeR$MacPolarizeR <- polarization_results

  return(object)
}

