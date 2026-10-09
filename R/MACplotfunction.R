#' Barplot of Macrophage Polarization
#'
#' @param object A Seurat object
#' @param group.by Metadata column for grouping (default: from CreateMacObj)
#' @param base.size Font size
#' @param cols Custom colors
#' @import ggplot2
#' @export
MacBarplot <- function(object, group.by = NULL, base.size = 14, cols = c("#C5283D", "#E9724C", "#FFC857")) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  if (is.null(object$MacPolarizeR)) {
    stop("MacPolarizeR prediction not found. Run 'MacPolarizeR' first.")
  }

  df <- as.data.frame(table(object@meta.data[[group.by]], object$MacPolarizeR))
  colnames(df) <- c("Group", "Polarization", "Count")

  # Calculate percentages
  df <- do.call(rbind, lapply(split(df, df$Group), function(x) {
    x$Percentage <- x$Count / sum(x$Count) * 100
    x
  }))

  p <- ggplot(df, aes(x = Group, y = Percentage, fill = Polarization)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Percentage", fill = "") +
    theme_classic(base_size = base.size) +
    scale_fill_manual(values = cols)

  return(p)
}

#' Radar Plot of Macrophage Polarization
#'
#' @param object A Seurat object
#' @param group.by Metadata column for grouping (default: from CreateMacObj)
#' @param split.by Facet column
#' @param sort Sort by count
#' @param base.size Font size
#' @param cols Custom colors
#' @param ncol Columns for facets
#' @import ggplot2
#' @export
MacRadar <- function(object, group.by = NULL, split.by = NULL, sort = FALSE, base.size = 14, cols = c("#C5283D", "#E9724C", "#FFC857"), ncol = 3) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  if (is.null(object$MacPolarizeR)) {
    stop("MacPolarizeR prediction not found. Run 'MacPolarizeR' first.")
  }

  dfplot <- object@meta.data
  if (sort) {
    dfplot[[group.by]] <- factor(dfplot[[group.by]], levels = rev(names(sort(table(dfplot[[group.by]])))), ordered = TRUE)
  }

  p <- ggplot(dfplot, aes(x = .data[[group.by]], fill = MacPolarizeR)) +
    geom_bar() +
    coord_polar() +
    labs(x = "", y = "", fill = "") +
    theme_void(base_size = base.size) +
    theme(axis.text.x = element_text(colour = "black")) +
    scale_fill_manual(values = cols)

  if (!is.null(split.by)) {
    p <- p + facet_wrap(vars(.data[[split.by]]), ncol = ncol)
  }

  return(p)
}

#' Polarization Cartesian Plot
#'
#' @param object A Seurat object
#' @param group.by Metadata column for coloring (default: from CreateMacObj)
#' @param split.by Facet column
#' @param style Plot style: "density" or "point"
#' @param fill Use fill for density
#' @param alpha Alpha for fill
#' @param cols Custom colors
#' @param base.size Font size
#' @param ncol Columns for facets
#' @import ggplot2
#' @export
PolCart <- function(object, group.by = NULL, split.by = NULL, style = "density", fill = TRUE, alpha = 0.1, cols = NULL, base.size = 15, ncol = 3) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (is.null(group.by)) {
    group.by <- object@misc$MACanalyzeR$ident
  }

  pol_data <- object@misc$MACanalyzeR$MacPolarizeR$SingleCell
  if (is.null(pol_data)) {
    stop("Run 'MacPolarizeR' first")
  }

  dfplot <- cbind(object@meta.data, pol_data)

  p <- ggplot(dfplot, aes(x = M1, y = M2, color = .data[[group.by]], fill = .data[[group.by]])) +
    labs(x = "M1 Score", y = "M2 Score", color = group.by, fill = group.by) +
    geom_vline(xintercept = 1, colour = 'grey30', linetype = "dashed") +
    geom_hline(yintercept = 1, colour = 'grey30', linetype = "dashed") +
    theme_classic(base_size = base.size)

  if (!is.null(cols)) {
    p <- p + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
  }

  if (style == 'point') {
    p <- p + geom_point(size = 1, alpha = 0.6)
  } else {
    p <- p + stat_density_2d(geom = "polygon", alpha = alpha, position = "identity")
  }

  if (!is.null(split.by)) {
    p <- p + facet_wrap(vars(.data[[split.by]]), ncol = ncol)
  }

  return(p)
}
