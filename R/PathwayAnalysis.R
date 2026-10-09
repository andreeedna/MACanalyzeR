#' Pathway Analysis for Macrophages
#'
#' @param object A Seurat object
#' @param pathways A list of pathways (optional)
#' @param group.by Metadata column to use for grouping (default: from CreateMacObj)
#' @param assay Assay to use (default: DefaultAssay(object))
#' @param layer Layer to use (default: "data")
#' @return A Seurat object with pathway scores
#' @importFrom stats weighted.mean var wilcox.test
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom future.apply future_lapply
#' @importFrom Seurat GetAssayData DefaultAssay
#' @export
PathAnalyzeR <- function(
  object,
  pathways = NULL,
  group.by = NULL,
  assay = NULL,
  layer = "data"
) {
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

  if (is.null(pathways)) {
    if (mac_info$organism == "mm") {
      pathways <- MACanalyzeR::mouse_selected
    } else {
      pathways <- MACanalyzeR::human_selected
    }
  } else if (!is.list(pathways)) {
    stop("Pathways must be a list of gene sets")
  }

  path_names <- names(pathways)
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  cond_tot <- object@meta.data[[group.by]]
  cond <- levels(factor(cond_tot))

  check.cols <- function(vect) {
    all(vect >= 0.001)
  }

  score_pathway <- function(p_name, store_mode) {
    genes <- pathways[[p_name]]
    genes <- intersect(genes, rownames(mtx))

    if (length(genes) < 5) {
      return(NULL)
    }

    mtx_pathway <- mtx[genes, , drop = FALSE]
    mtx_pathway <- mtx_pathway[Matrix::rowSums(mtx_pathway) > 0, , drop = FALSE]

    # Use a dense representation for the small subset of genes
    mtx_dense <- as.matrix(mtx_pathway)
    sample_mean <- apply(mtx_dense, 1, function(x) by(x, cond_tot, mean))

    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if (length(keep) < 3) {
      return(NULL)
    }

    mtx_dense <- mtx_dense[keep, , drop = FALSE]
    mean_expr_by_cluster <- apply(mtx_dense, 1, function(x) {
      by(x, cond_tot, mean)
    })

    if (store_mode == "SingleCell") {
      ratio_expr <- mtx_dense / colMeans(mean_expr_by_cluster)
      gene_weight <- apply(mtx_dense, 1, stats::var)
      mean_exp_pathway <- apply(ratio_expr, 2, function(x) {
        stats::weighted.mean(x, gene_weight / sum(gene_weight))
      })
      return(mean_exp_pathway)
    } else {
      ratio_expr <- t(mean_expr_by_cluster) / colMeans(mean_expr_by_cluster)
      gene_weight <- apply(mtx_dense, 1, stats::var)
      mean_exp_pathway <- apply(ratio_expr, 2, function(x) {
        stats::weighted.mean(x, gene_weight / sum(gene_weight))
      })
      pval <- stats::wilcox.test(ratio_expr[, 1], mu = 1)$p.value # Note: original logic only tests ratio_expr[,1]
      return(c(mean_exp_pathway[cond], pval = pval))
    }
  }

  results <- list()
  for (store in c("Total", "SingleCell")) {
    message("Calculating ", store, " scores...")
    # Parallel processing using future_lapply
    scored_list <- future.apply::future_lapply(
      path_names,
      score_pathway,
      store_mode = store,
      future.seed = TRUE
    )
    names(scored_list) <- path_names

    # Filter out NULLs
    scored_list <- scored_list[!sapply(scored_list, is.null)]

    if (store == "SingleCell") {
      res_df <- as.data.frame(do.call(cbind, scored_list))
      # Add to meta.data
      object@meta.data <- cbind(object@meta.data, res_df)
    } else {
      res_df <- as.data.frame(do.call(rbind, scored_list))
      colnames(res_df) <- c(cond, "pval")
    }
    results[[store]] <- res_df
  }

  object@misc$MACanalyzeR$PathAnalyzeR[[group.by]] <- results
  return(object)
}

#' Heatmap of Pathway Scores
#'
#' @param object A Seurat object
#' @param group.by Metadata column used for PathAnalyzeR
#' @param cols Colors for heatmap
#' @param pval P-value cutoff
#' @param base.size Font size
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
PathHeat <- function(
  object,
  group.by = NULL,
  cols = c("blue", "white", "red"),
  pval = 0.05,
  base.size = 10
) {
  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  dfplot <- object@misc$MACanalyzeR$PathAnalyzeR[[group.by]][["Total"]]
  if (is.null(dfplot)) {
    stop("Run 'PathAnalyzeR' first for ", group.by)
  }

  dfplot <- dfplot[
    which(dfplot$pval < pval),
    colnames(dfplot) != "pval",
    drop = FALSE
  ]

  sort_row <- c()
  for (i in colnames(dfplot)) {
    select_row <- which(apply(dfplot, 1, max) == dfplot[, i])
    tmp <- rownames(dfplot)[select_row][order(
      dfplot[select_row, i],
      decreasing = TRUE
    )]
    sort_row <- c(sort_row, tmp)
  }
  dfplot[is.na(dfplot)] <- 1
  b <- max((max(dfplot) - 1), (1 - min(dfplot)))
  mybreaks <- c(
    seq((1 - b), 1 - (b / 3), length.out = 33),
    seq(1 - (b / 3) + 0.01, 1 + (b / 3), length.out = 34),
    seq(1 + (b / 3) + 0.01, 1 + b, length.out = 33)
  )
  color <- grDevices::colorRampPalette(cols)(100)
  pheatmap::pheatmap(
    dfplot[sort_row, ],
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = color,
    breaks = mybreaks,
    fontsize = base.size,
    border_color = FALSE
  )
}

#' Display Pathways Used
#'
#' @param object A Seurat object
#' @param group.by Metadata column used
#' @export
PathDisplay <- function(object, group.by = NULL) {
  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }
  paths <- colnames(object@misc$MACanalyzeR$PathAnalyzeR[[group.by]][[
    "SingleCell"
  ]])
  if (is.null(paths)) {
    stop(
      "No pathway results found for ",
      group.by,
      ". Run 'PathAnalyzeR' first."
    )
  }
  message("Pathways used for ", group.by, ":")
  for (i in seq_along(paths)) {
    message(i, "-", paths[i])
  }
}

#' Cartesian Plot of Pathways
#'
#' @param object A Seurat object
#' @param pathways Indices or names of pathways to plot (1 or 2)
#' @param group.by Metadata column for coloring
#' @param split.by Metadata column for facetting
#' @param fill Use fill for density
#' @param alpha Alpha for fill
#' @param cols Custom colors
#' @param base.size Font size
#' @param ncol Number of columns for facetting
#' @import ggplot2
#' @export
PathCart <- function(
  object,
  pathways,
  group.by = NULL,
  split.by = NULL,
  fill = TRUE,
  alpha = 0.1,
  cols = NULL,
  base.size = 15,
  ncol = 3
) {
  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  all_pathway <- colnames(object@misc$MACanalyzeR$PathAnalyzeR[[group.by]][[
    "SingleCell"
  ]])
  if (is.null(all_pathway)) {
    stop("Run 'PathAnalyzeR' first")
  }

  if (is.numeric(pathways)) {
    if (max(pathways) > length(all_pathway)) {
      stop("Pathway index out of range")
    }
    path_names <- all_pathway[pathways]
  } else {
    if (!all(pathways %in% all_pathway)) {
      stop("One or more pathways not found in results")
    }
    path_names <- pathways
  }

  dfplot <- cbind(
    object@meta.data[, path_names, drop = FALSE],
    object@meta.data
  )

  if (length(path_names) == 1) {
    plot <- ggplot(
      dfplot,
      aes(x = .data[[path_names[1]]], color = .data[[group.by]])
    ) +
      labs(x = path_names[1], y = "", color = group.by, fill = group.by) +
      geom_vline(xintercept = 1, colour = 'grey30') +
      theme_classic(base_size = base.size)

    if (fill) {
      plot <- plot + geom_density(aes(fill = .data[[group.by]]), alpha = alpha)
    } else {
      plot <- plot + geom_density()
    }
  } else if (length(path_names) == 2) {
    plot <- ggplot(
      dfplot,
      aes(x = .data[[path_names[1]]], y = .data[[path_names[2]]])
    ) +
      labs(
        x = path_names[1],
        y = path_names[2],
        color = group.by,
        fill = group.by
      ) +
      geom_vline(xintercept = 1, colour = 'grey30') +
      geom_hline(yintercept = 1, colour = 'grey30') +
      theme_classic(base_size = base.size)

    if (fill) {
      plot <- plot +
        stat_density_2d(
          aes(color = .data[[group.by]], fill = .data[[group.by]]),
          geom = "polygon",
          alpha = alpha,
          position = "identity"
        )
    } else {
      plot <- plot + geom_density_2d(aes(color = .data[[group.by]]))
    }
  } else {
    stop("Please specify one or two pathways")
  }

  if (!is.null(cols)) {
    plot <- plot +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols)
  }

  if (!is.null(split.by)) {
    plot <- plot + facet_wrap(vars(.data[[split.by]]), ncol = ncol)
  }

  return(plot)
}
