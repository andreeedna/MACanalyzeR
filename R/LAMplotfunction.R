.check_foam <- function(mac_obj) {
  if (is.null(mac_obj@MetaData[["fMAC+"]])) {
    stop("FoamSpotteR prediction is absent. Run FoamSpotteR() first.", call. = FALSE)
  }
}

FoamPlot <- function(mac_obj, txt.size = 14, shade = "viridis", split.by = NULL, pt.size = 2, ncol = 3) {
  .check_obj(mac_obj)
  .check_foam(mac_obj)
  .check_col(mac_obj, split.by, "split.by")

  emb <- .embedding(mac_obj)
  dfplot <- cbind(mac_obj@MetaData, emb$data)
  dfplot <- dfplot[order(dfplot[["fMAC+"]], decreasing = FALSE), , drop = FALSE]

  plot <- ggplot(dfplot, aes(x = .data[[".dim1"]], y = .data[[".dim2"]], color = .data[["fMAC+"]])) +
    labs(x = emb$labels[1], y = emb$labels[2], color = "FoamDEX") +
    geom_point(size = pt.size) +
    .mac_theme(txt.size) +
    scale_color_viridis_c(option = shade, limits = c(0, 1))

  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  plot
}


###################################

FoamLine <- function(mac_obj, plot.by = mac_obj@ident, split.by = NULL, txt.size = 14, fill = TRUE, a = 0.5,
                     col = NULL, ncol = 3) {
  .check_obj(mac_obj)
  .check_foam(mac_obj)
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")

  dfplot <- mac_obj@MetaData
  plot <- ggplot(dfplot, aes(x = .data[["fMAC+"]], color = .data[[plot.by]])) +
    labs(x = "FoamDEX", y = "", fill = plot.by, color = plot.by) +
    xlim(0, 1) +
    .mac_theme(txt.size)

  if (fill) {
    plot <- plot + geom_density(aes(y = after_stat(scaled), fill = .data[[plot.by]]), alpha = a)
  } else {
    plot <- plot + geom_density(aes(y = after_stat(scaled)))
  }

  if (.check_colors(col, dfplot[[plot.by]], "plot.by")) {
    plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
  }
  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  plot
}
