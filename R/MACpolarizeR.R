## Polarization classes, from the most M2-like to the most M1-like.
.pol_levels <- c("Inflammatory", "Transitional", "Healing")

## Cluster the polarization balance into 3 states and name them by centroid:
## lowest balance -> Healing, highest -> Inflammatory, middle -> Transitional.
.pol_classify <- function(balance, nstart) {
  km <- stats::kmeans(balance, centers = 3, nstart = nstart, iter.max = 100)
  rank <- rank(km$centers[, 1], ties.method = "first")
  factor(c("Healing", "Transitional", "Inflammatory")[rank[km$cluster]], levels = .pol_levels)
}

## Per-cell M1 and M2 scores (MACanalyzeR score relative to the whole dataset)
## and their standardized difference (polarization balance).
.pol_scores <- function(m1, m2) {
  one <- rep("all", ncol(m1))
  s1 <- .score_set(m1, one)
  s2 <- .score_set(m2, one)
  if (is.null(s1) || is.null(s2)) return(NULL)
  sc <- cbind(M1 = s1$single, M2 = s2$single)
  z <- scale(sc)
  list(scores = sc, balance = as.vector(z[, "M1"] - z[, "M2"]), genes = list(M1 = s1$genes, M2 = s2$genes))
}


MacPolarizeR <- function(mac_obj, seed = 42, nstart = 25, n.boot = 20) {
  .check_obj(mac_obj)
  inf <- .org_data(mac_obj, MACanalyzeR::mouse_PolGenes, MACanalyzeR::human_PolGenes)

  m1 <- .expr(mac_obj, inf$M1)
  m2 <- .expr(mac_obj, inf$M2)
  pol <- .pol_scores(m1, m2)
  if (is.null(pol)) {
    stop("Too few M1/M2 marker genes are expressed (found ", nrow(m1), " M1 and ", nrow(m2),
         " M2 genes; at least 5 of each, with 3 expressed on average, are needed)", call. = FALSE)
  }
  if (stats::sd(pol$balance) == 0) stop("The polarization balance is constant across cells", call. = FALSE)

  res <- .with_seed(seed, {
    mac <- .pol_classify(pol$balance, nstart)

    ## stability: resample the marker genes and repeat the classification
    stab <- NULL
    if (n.boot > 0) {
      m1u <- m1[pol$genes$M1, , drop = FALSE]
      m2u <- m2[pol$genes$M2, , drop = FALSE]
      agree <- vapply(seq_len(n.boot), function(b) {
        pb <- .pol_scores(m1u[sample(nrow(m1u), replace = TRUE), , drop = FALSE],
                          m2u[sample(nrow(m2u), replace = TRUE), , drop = FALSE])
        if (is.null(pb) || stats::sd(pb$balance) == 0) return(rep(NA, length(mac)))
        as.character(.pol_classify(pb$balance, nstart)) == as.character(mac)
      }, logical(length(mac)))
      stab <- rowMeans(matrix(agree, nrow = length(mac)), na.rm = TRUE)
    }
    list(mac = mac, stab = stab)
  })

  total <- t(apply(pol$scores, 2, function(x) tapply(x, res$mac, mean)))
  total <- total[, .pol_levels, drop = FALSE]

  mac_obj@MetaData$Mac <- res$mac
  if (!is.null(res$stab)) mac_obj@MetaData$Mac.stability <- res$stab

  mac_obj@MacPolarizeR <- list(SingleCell = pol$scores,
                               Total = total,
                               Balance = stats::setNames(pol$balance, colnames(m1)),
                               Genes = pol$genes,
                               Params = list(seed = seed, nstart = nstart, n.boot = n.boot,
                                             version = as.character(utils::packageVersion("MACanalyzeR"))))

  msg <- paste0("MacPolarizeR: ", paste(names(table(res$mac)), table(res$mac), sep = " = ", collapse = ", "))
  if (!is.null(res$stab)) msg <- paste0(msg, " | mean stability ", round(mean(res$stab, na.rm = TRUE), 2))
  message(msg)
  mac_obj
}


.pol_colors <- c("Inflammatory" = "#C5283D", "Transitional" = "#E9724C", "Healing" = "#FFC857")

.pol_fill <- function(col) {
  if (length(col) == 3) {
    if (is.null(names(col))) names(col) <- .pol_levels
    return(col)
  }
  message("The colour vector must contain 3 colours: using the default palette")
  .pol_colors
}


MacBarplot <- function(mac_obj, plot.by = mac_obj@ident, txt.size = 14, col = c("#C5283D", "#E9724C", "#FFC857")) {
  .check_obj(mac_obj)
  .check_col(mac_obj, "Mac", "Mac")
  .check_col(mac_obj, plot.by, "plot.by")

  count <- table(mac_obj@MetaData[[plot.by]], factor(mac_obj@MetaData$Mac, levels = .pol_levels))
  perc <- prop.table(count, 1) * 100
  dfplot <- as.data.frame(perc, stringsAsFactors = FALSE)
  colnames(dfplot) <- c("group", "Mac", "value")
  dfplot$group <- factor(dfplot$group, levels = rownames(perc))
  dfplot$Mac <- factor(dfplot$Mac, levels = .pol_levels)

  ggplot(dfplot, aes(x = .data[["group"]], y = .data[["value"]], fill = .data[["Mac"]])) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "% of cells", fill = "") +
    .mac_theme(txt.size) +
    theme(axis.text.x = element_text(color = "black", size = txt.size),
          axis.text.y = element_text(color = "black", size = txt.size)) +
    scale_fill_manual(values = .pol_fill(col))
}



MacRadar <- function(mac_obj, plot.by = mac_obj@ident, split.by = NULL, sort = FALSE, txt.size = 14,
                     col = c("#C5283D", "#E9724C", "#FFC857"), ncol = 3) {
  .check_obj(mac_obj)
  .check_col(mac_obj, "Mac", "Mac")
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")

  dfplot <- mac_obj@MetaData
  tab <- table(dfplot[[plot.by]])
  if (sort) {
    dfplot[[plot.by]] <- factor(dfplot[[plot.by]], levels = names(sort(tab, decreasing = TRUE)))
  } else {
    dfplot[[plot.by]] <- factor(dfplot[[plot.by]])
  }
  n_lev <- nlevels(dfplot[[plot.by]])

  plot <- ggplot(dfplot, aes(x = .data[[plot.by]], fill = .data[["Mac"]])) +
    geom_bar() +
    coord_polar(start = -(pi / n_lev)) +
    labs(fill = "") +
    ylim(c(-max(tab) / 5, max(tab))) +
    theme_void() +
    theme(text = element_text(size = txt.size),
          axis.text.x = element_text(colour = "black")) +
    scale_fill_manual(values = .pol_fill(col))

  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  plot
}


PolCart <- function(mac_obj, plot.by = mac_obj@ident, split.by = NULL, style = "density", fill = TRUE,
                    a = 0.1, col = NULL, txt.size = 15, ncol = 3, scale = TRUE) {
  .check_obj(mac_obj)
  if (is.null(mac_obj@MacPolarizeR$SingleCell)) stop("Run MacPolarizeR() first.", call. = FALSE)
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")
  if (!style %in% c("density", "point")) stop("'style' must be 'density' or 'point'", call. = FALSE)

  sc <- mac_obj@MacPolarizeR$SingleCell
  dfplot <- mac_obj@MetaData
  dfplot$.M1 <- sc[, "M1"]
  dfplot$.M2 <- sc[, "M2"]
  lg <- .log_scores(dfplot, c(".M1", ".M2"), scale)
  dfplot <- lg$data

  plot <- ggplot(dfplot, aes(x = .data[[".M1"]], y = .data[[".M2"]],
                             color = .data[[plot.by]], fill = .data[[plot.by]])) +
    labs(x = .axis_lab("M1 score", scale), y = .axis_lab("M2 score", scale), color = plot.by, fill = plot.by) +
    geom_vline(xintercept = lg$ref, colour = "grey30") +
    geom_hline(yintercept = lg$ref, colour = "grey30") +
    .mac_theme(txt.size)

  if (style == "point") {
    plot <- plot + geom_point(size = 2)
  } else if (fill) {
    plot <- plot + stat_density_2d(geom = "polygon", alpha = a, position = "identity")
  } else {
    plot <- plot + geom_density_2d()
  }

  if (.check_colors(col, dfplot[[plot.by]], "plot.by")) {
    plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
  }
  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)
  plot
}
