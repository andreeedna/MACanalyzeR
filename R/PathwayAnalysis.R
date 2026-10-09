PathAnalyzeR <- function(mac_obj, pathway = NULL, meta = NULL, verbose = TRUE) {
  .check_obj(mac_obj)
  if (is.null(meta)) meta <- mac_obj@ident
  .check_col(mac_obj, meta, "meta")

  if (is.null(pathway)) {
    pathway <- .org_data(mac_obj, MACanalyzeR::mouse_selected, MACanalyzeR::human_selected)
  } else if (is.character(pathway)) {
    pathway <- list("Custom gene set" = pathway)
  } else if (!is.list(pathway)) {
    stop("'pathway' must be a named list of gene vectors. See MACanalyzeR::mouse_selected for an example.",
         call. = FALSE)
  }

  groups <- mac_obj@MetaData[[meta]]
  if (length(unique(stats::na.omit(groups))) < 2) {
    warning("'", meta, "' has only one group: group-level scores are all 1 and no p-value is computed",
            call. = FALSE)
  }

  res <- .score_list(mac_obj, pathway, groups, single = TRUE, test = TRUE, verbose = verbose,
                     label = paste0("PathAnalyzeR on ", meta, " (", length(pathway), " gene sets)"))

  total <- as.data.frame(res$total)
  total$pval <- res$pval
  total$padj <- stats::p.adjust(res$pval, method = "BH")

  mac_obj@PathAnalyzeR[[meta]] <- list(Total = total, SingleCell = as.data.frame(res$single, check.names = FALSE))
  mac_obj
}

#####

PathHeat <- function(mac_obj, plot.by = mac_obj@ident, col = c("blue", "white", "red"), pval = 0.05,
                     txt.size = 10, adjust = FALSE) {
  .check_obj(mac_obj)
  dfplot <- mac_obj@PathAnalyzeR[[plot.by]][["Total"]]
  if (is.null(dfplot)) {
    stop("No PathAnalyzeR results for '", plot.by, "'. Run PathAnalyzeR(mac_obj, meta = '", plot.by, "') first.",
         call. = FALSE)
  }

  pcol <- if (adjust && "padj" %in% colnames(dfplot)) "padj" else "pval"
  groups <- setdiff(colnames(dfplot), c("pval", "padj"))
  if (all(is.na(dfplot[[pcol]]))) {
    keep <- rowSums(!is.na(dfplot[, groups, drop = FALSE])) > 0
    message("No p-values available: showing all scored pathways")
  } else {
    keep <- !is.na(dfplot[[pcol]]) & dfplot[[pcol]] < pval
  }
  dfplot <- as.matrix(dfplot[keep, groups, drop = FALSE])
  if (nrow(dfplot) == 0) {
    message("No pathway with ", pcol, " < ", pval, ": nothing to plot")
    return(invisible(NULL))
  }

  dfplot <- dfplot[.heat_order(dfplot), , drop = FALSE]
  dfplot[is.na(dfplot)] <- 1
  pheatmap(dfplot, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(col)(100),
           breaks = .heat_breaks(dfplot), fontsize = txt.size, border_color = NA)
}


PathDisplay <- function(mac_obj, meta = mac_obj@ident) {
  .check_obj(mac_obj)
  res <- .single_cell_results(mac_obj, meta)
  message("Pathways used for the ", meta, " pathway analysis:")
  p <- colnames(res)
  message(paste0(seq_along(p), " - ", p, collapse = "\n"))
  invisible(data.frame(index = seq_along(p), pathway = p))
}


PathCart <- function(mac_obj, pathway, plot.by = mac_obj@ident, split.by = NULL, fill = TRUE, a = 0.1,
                     col = NULL, txt.size = 15, ncol = 3, scale = TRUE) {
  .check_obj(mac_obj)
  res <- .single_cell_results(mac_obj, plot.by)
  .check_col(mac_obj, split.by, "split.by")
  pathway <- .pick_pathways(colnames(res), pathway)
  if (length(pathway) > 2) stop("Please select one or two pathways", call. = FALSE)

  dfplot <- mac_obj@MetaData
  dfplot$.p1 <- res[[pathway[1]]]
  if (length(pathway) == 2) dfplot$.p2 <- res[[pathway[2]]]
  ok <- !is.na(dfplot$.p1)
  if (length(pathway) == 2) ok <- ok & !is.na(dfplot$.p2)
  dfplot <- dfplot[ok, , drop = FALSE]
  if (nrow(dfplot) == 0) stop("The selected pathway(s) could not be scored", call. = FALSE)

  lg <- .log_scores(dfplot, c(".p1", ".p2")[seq_along(pathway)], scale)
  dfplot <- lg$data
  if (nrow(dfplot) == 0) stop("No cell left to plot", call. = FALSE)

  if (length(pathway) == 1) {
    plot <- ggplot(dfplot, aes(x = .data[[".p1"]], color = .data[[plot.by]])) +
      labs(x = .axis_lab(pathway[1], scale), y = "", color = "", fill = "") +
      geom_vline(xintercept = lg$ref, colour = "grey30")
    plot <- if (fill) plot + geom_density(aes(fill = .data[[plot.by]]), alpha = a) else plot + geom_density()
  } else {
    plot <- ggplot(dfplot, aes(x = .data[[".p1"]], y = .data[[".p2"]])) +
      labs(x = .axis_lab(pathway[1], scale), y = .axis_lab(pathway[2], scale),
           color = plot.by, fill = plot.by) +
      geom_vline(xintercept = lg$ref, colour = "grey30") +
      geom_hline(yintercept = lg$ref, colour = "grey30")
    plot <- if (fill) {
      plot + stat_density_2d(aes(color = .data[[plot.by]], fill = .data[[plot.by]]),
                             geom = "polygon", alpha = a, position = "identity")
    } else {
      plot + geom_density_2d(aes(color = .data[[plot.by]]))
    }
  }

  plot <- plot + .mac_theme(txt.size)
  if (.check_colors(col, dfplot[[plot.by]], "plot.by")) {
    plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
  }
  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  plot
}


PathPlot <- function(mac_obj, pathway, txt.size = 14, shade = "turbo", split.by = NULL, pt.size = 2,
                     ncol = 2, max.cutoff = NA, meta = mac_obj@ident) {
  .check_obj(mac_obj)
  res <- .single_cell_results(mac_obj, meta)
  .check_col(mac_obj, split.by, "split.by")
  pathway <- .pick_pathways(colnames(res), pathway)

  emb <- .embedding(mac_obj)
  base <- cbind(mac_obj@MetaData, emb$data)

  ## long format: one block of cells per pathway
  dfplot <- do.call(rbind, lapply(pathway, function(p) {
    d <- base
    d$.score <- res[[p]]
    d$.pathway <- p
    d
  }))
  dfplot$.pathway <- factor(dfplot$.pathway, levels = pathway)

  if (!is.na(max.cutoff)) {
    if (max.cutoff < 2) warning("A cutoff below 2 is not recommended", call. = FALSE)
    dfplot$.score <- pmin(dfplot$.score, max.cutoff)
  }
  dfplot <- dfplot[order(dfplot$.score, na.last = FALSE), , drop = FALSE]

  plot <- ggplot(dfplot, aes(x = .data[[".dim1"]], y = .data[[".dim2"]], color = .data[[".score"]])) +
    geom_point(size = pt.size) +
    labs(x = emb$labels[1], y = emb$labels[2], color = "Score") +
    .mac_theme(txt.size) +
    scale_color_viridis_c(option = shade)

  if (length(pathway) == 1) {
    plot <- plot + ggtitle(pathway)
    if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  } else if (is.null(split.by)) {
    plot <- plot + facet_wrap(".pathway", ncol = ncol)
  } else {
    plot <- plot + facet_grid(stats::as.formula(paste0(".pathway ~ `", split.by, "`")))
  }
  plot
}

###############################################

PathViolin <- function(mac_obj, pathway, plot.by = mac_obj@ident, txt.size = 15, ncol = 1, col = NULL,
                       intercept = FALSE, stat = FALSE) {
  .check_obj(mac_obj)
  res <- .single_cell_results(mac_obj, plot.by)
  pathway <- .pick_pathways(colnames(res), pathway)
  if (is.null(intercept)) intercept <- FALSE

  dfplot <- mac_obj@MetaData
  dfplot[[plot.by]] <- factor(dfplot[[plot.by]])
  use_col <- .check_colors(col, dfplot[[plot.by]], "plot.by")

  plot <- lapply(pathway, function(p) {
    d <- dfplot
    d$.score <- res[[p]]
    d <- d[!is.na(d$.score) & !is.na(d[[plot.by]]), , drop = FALSE]

    p1 <- ggplot(d, aes(x = .data[[plot.by]], y = .data[[".score"]],
                        fill = .data[[plot.by]], color = .data[[plot.by]])) +
      geom_violin() +
      geom_boxplot(width = 0.3, color = "grey30", alpha = 0.2) +
      labs(x = "", y = p, color = plot.by, fill = plot.by) +
      .mac_theme(txt.size) +
      theme(legend.position = "none")

    if (isTRUE(intercept)) p1 <- p1 + geom_hline(yintercept = 1, colour = "grey30")
    if (use_col) p1 <- p1 + scale_fill_manual(values = col) + scale_color_manual(values = col)

    if (stat) {
      lev <- levels(droplevels(d[[plot.by]]))
      if (length(lev) > 1) {
        comp <- utils::combn(lev, 2, simplify = FALSE)
        p1 <- p1 + geom_signif(comparisons = comp, map_signif_level = FALSE, color = "black",
                               step_increase = 0.1)
      }
    }
    p1
  })
  names(plot) <- pathway

  grid.arrange(grobs = plot, ncol = ncol)
}
