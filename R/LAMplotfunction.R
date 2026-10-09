#' Foam Cell Density Plot
#'
#' @param object A Seurat object
#' @param group.by Metadata column for coloring (default: from CreateMacObj)
#' @param split.by Metadata column for facetting
#' @param base.size Font size
#' @param fill Use fill for density
#' @param alpha Alpha for fill
#' @param cols Custom colors
#' @param ncol Columns for facetting
#' @import ggplot2
#' @export
FoamLine <- function(object, group.by = NULL, split.by = NULL, base.size = 14, fill = TRUE, alpha = 0.5, cols = NULL, ncol = 3) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  if (is.null(object$Foam_fMAC.)) {
    stop("FoamSpotteR prediction not found. Run 'FoamSpotteR' first.")
  }

  dfplot <- object@meta.data

  p <- ggplot(dfplot, aes(x = Foam_fMAC., color = .data[[group.by]])) +
    labs(x = "FoamDEX", y = 'Density', fill = group.by, color = group.by) +
    xlim(0, 1) +
    theme_classic(base_size = base.size)

  if (fill) {
    p <- p + geom_density(aes(fill = .data[[group.by]]), alpha = alpha)
  } else {
    p <- p + geom_density()
  }

  if (!is.null(cols)) {
    p <- p + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
  }

  if (!is.null(split.by)) {
    p <- p + facet_wrap(vars(.data[[split.by]]), ncol = ncol)
  }

  return(p)
}



