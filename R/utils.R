## Internal helpers shared by the MACanalyzeR modules.
## Names start with a dot, so they are not exported by exportPattern().

.check_obj <- function(mac_obj) {
  if (!methods::is(mac_obj, "MACanalyzeR")) {
    stop("The input object is not a MACanalyzeR object. ",
         "Create it with the CreateMacObj() function.", call. = FALSE)
  }
  invisible(TRUE)
}

## Check that a metadata column exists, with a helpful message when it is
## produced by a module that has not been run yet.
.check_col <- function(mac_obj, col, arg = "plot.by") {
  if (is.null(col)) return(invisible(TRUE))
  if (length(col) != 1 || !is.character(col)) {
    stop("'", arg, "' must be the name of one column of mac_obj@MetaData", call. = FALSE)
  }
  if (!col %in% colnames(mac_obj@MetaData)) {
    hint <- switch(col,
                   Mac  = "Run MacPolarizeR() first.",
                   Foam = "Run FoamSpotteR() first.",
                   paste0("Available columns: ", paste(colnames(mac_obj@MetaData), collapse = ", ")))
    where <- if (identical(arg, col)) "" else paste0(" (argument '", arg, "')")
    stop("Column '", col, "'", where, " is not in mac_obj@MetaData. ", hint, call. = FALSE)
  }
  invisible(TRUE)
}

## Check that a colour vector has one colour per group.
.check_colors <- function(col, groups, arg = "plot.by") {
  if (is.null(col)) return(FALSE)
  n <- length(unique(stats::na.omit(groups)))
  if (length(col) != n) {
    stop("'col' must contain one colour per group of ", arg, " (", n, " groups), not ", length(col), call. = FALSE)
  }
  TRUE
}

## Organism-specific reference data.
.org_data <- function(mac_obj, mouse, human) {
  if (identical(mac_obj@organism, "hs")) human else mouse
}

## Dense expression matrix (genes x cells) for a set of genes.
## Works both with a dense matrix and with a sparse Matrix object.
.expr <- function(mac_obj, genes) {
  m <- mac_obj@CountData@NormalizedCount
  genes <- intersect(unique(genes), rownames(m))
  as.matrix(m[genes, , drop = FALSE])
}

## Group means of a genes x cells matrix: genes x groups.
## Cells with NA group are ignored (as tapply/by do).
.group_means <- function(m, groups) {
  groups <- factor(groups)
  lev <- levels(groups)
  ind <- matrix(0, nrow = length(groups), ncol = length(lev), dimnames = list(NULL, lev))
  ok <- !is.na(groups)
  ind[cbind(which(ok), as.integer(groups[ok]))] <- 1
  sweep(m %*% ind, 2, colSums(ind), "/")
}

.row_vars <- function(m) {
  n <- ncol(m)
  mu <- rowMeans(m)
  (rowSums(m^2) - n * mu^2) / (n - 1)
}

## Core MACanalyzeR gene set score.
##
## For a gene set, the expression of each gene is divided by its average
## expression across groups (mean of the group means); ratios are then averaged
## across genes, weighting each gene by its expression variance across cells.
## A score of 1 means average activity.
##
## Gene filters (as in the original implementation): at least `min.genes` genes
## of the set must be detected in the dataset; genes never expressed are dropped;
## genes with a mean below `min.mean` in at least one group are dropped and at
## least `min.keep` genes must remain. Otherwise NULL is returned.
##
## Returns a list with:
##   total  - named vector, one score per group
##   single - one score per cell (NULL if single = FALSE)
##   ratio  - genes x groups matrix of gene ratios (used for testing)
##   genes  - genes used
.score_set <- function(m, groups, single = TRUE, min.genes = 5, min.keep = 3, min.mean = 0.001) {
  if (nrow(m) < min.genes) return(NULL)
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  gm <- .group_means(m, groups)
  keep <- apply(gm, 1, function(x) all(x >= min.mean))
  if (sum(keep) < min.keep) return(NULL)
  m <- m[keep, , drop = FALSE]
  gm <- gm[keep, , drop = FALSE]

  ref <- rowMeans(gm)
  w <- .row_vars(m)
  w <- w / sum(w)

  ratio <- gm / ref
  total <- colSums(ratio * w)

  sc <- NULL
  if (single) sc <- colSums((m / ref) * w)

  list(total = total, single = sc, ratio = ratio, genes = rownames(m))
}

## P-value for differences between groups of a gene set score.
## Two groups: one-sample Wilcoxon signed-rank test of the gene ratios of the
## first group against 1 (with two groups this is equivalent to testing the
## difference between them). Three or more groups: Friedman test with genes as
## blocks. One group: NA.
.score_pval <- function(ratio) {
  k <- ncol(ratio)
  if (k < 2) return(NA_real_)
  p <- tryCatch(suppressWarnings(
    if (k == 2) stats::wilcox.test(ratio[, 1], mu = 1)$p.value
    else stats::friedman.test(ratio)$p.value
  ), error = function(e) NA_real_)
  p
}

## Score a list of gene sets. Returns list(Total = data.frame, SingleCell = data.frame).
.score_list <- function(mac_obj, pathway, groups, single = TRUE, test = TRUE, verbose = TRUE, label = NULL) {
  path_names <- names(pathway)
  if (is.null(path_names) || any(path_names == "") || anyDuplicated(path_names)) {
    stop("'pathway' must be a list with unique, non-empty names", call. = FALSE)
  }
  lev <- levels(factor(groups))
  cells <- colnames(mac_obj@CountData@NormalizedCount)

  total <- matrix(NA_real_, nrow = length(path_names), ncol = length(lev),
                  dimnames = list(path_names, lev))
  pval <- stats::setNames(rep(NA_real_, length(path_names)), path_names)
  sc <- if (single) matrix(NA_real_, nrow = length(cells), ncol = length(path_names),
                           dimnames = list(cells, path_names)) else NULL

  if (verbose && !is.null(label)) message(label)
  pb <- NULL
  if (verbose && length(path_names) > 1) {
    pb <- utils::txtProgressBar(min = 0, max = length(path_names), style = 3, width = 60, char = "+")
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_along(path_names)) {
    p <- path_names[i]
    s <- .score_set(.expr(mac_obj, pathway[[p]]), groups, single = single)
    if (!is.null(s)) {
      total[p, ] <- s$total[lev]
      if (test) pval[p] <- .score_pval(s$ratio)
      if (single) sc[, p] <- s$single
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
  }

  n_na <- sum(is.na(total[, 1]))
  if (verbose && n_na > 0) {
    message(n_na, " of ", length(path_names), " gene sets were skipped (too few expressed genes)")
  }

  list(total = total, pval = pval, single = sc)
}

## Row order for heatmaps: pathways grouped by the column where they peak,
## then by decreasing value.
.heat_order <- function(dat) {
  dat <- as.matrix(dat)
  best <- apply(dat, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
  peak <- apply(dat, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  rownames(dat)[order(best, -peak, na.last = TRUE)]
}

## Colour breaks symmetric around 1 (average activity).
.heat_breaks <- function(dat) {
  b <- max(abs(range(dat, na.rm = TRUE) - 1))
  if (!is.finite(b) || b == 0) b <- 0.01
  seq(1 - b, 1 + b, length.out = 101)
}

## Run an expression with a fixed seed without changing the user's RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had_seed) get(".Random.seed", envir = genv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)
  expr
}

## Shared ggplot2 theme of the package.
.mac_theme <- function(txt.size) {
  theme(text = element_text(size = txt.size),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black"))
}

## Embedding columns of the object (first two dimensions).
.embedding <- function(mac_obj) {
  emb <- mac_obj@Reduction
  if (ncol(emb) < 2) stop("The reduction stored in the object has less than two dimensions", call. = FALSE)
  emb <- emb[, 1:2, drop = FALSE]
  colnames(emb) <- c(".dim1", ".dim2")
  list(data = emb, labels = colnames(mac_obj@Reduction)[1:2])
}

## Pathway indices or names -> pathway names.
.pick_pathways <- function(available, pathway) {
  if (is.character(pathway)) {
    miss <- setdiff(pathway, available)
    if (length(miss) > 0) stop("Pathway(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    return(pathway)
  }
  if (!is.numeric(pathway) || any(pathway < 1) || any(pathway > length(available))) {
    stop("'pathway' must contain indices between 1 and ", length(available),
         " (see PathDisplay()) or pathway names", call. = FALSE)
  }
  available[pathway]
}

.single_cell_results <- function(mac_obj, meta, fun = "PathAnalyzeR") {
  res <- mac_obj@PathAnalyzeR[[meta]][["SingleCell"]]
  if (is.null(res)) {
    stop("No PathAnalyzeR results for '", meta, "'. Run ", fun, "(mac_obj, meta = '", meta, "') first.", call. = FALSE)
  }
  res
}

utils::globalVariables("scaled")

## Optional log2 transformation of score columns for the Cartesian plots.
## Scores are ratios to the average (1 = average activity), so log2 makes
## increases and decreases symmetric around 0 (+1 = twice, -1 = half the average).
## Cells with a score <= 0 cannot be log-transformed and are removed.
.log_scores <- function(dfplot, cols, scale) {
  if (!is.logical(scale) || length(scale) != 1 || is.na(scale)) {
    stop("'scale' must be TRUE or FALSE", call. = FALSE)
  }
  if (!scale) return(list(data = dfplot, ref = 1))
  ok <- Reduce(`&`, lapply(cols, function(v) !is.na(dfplot[[v]]) & dfplot[[v]] > 0))
  if (any(!ok)) message("Removed ", sum(!ok), " cells with a score <= 0, which cannot be shown on a log scale")
  dfplot <- dfplot[ok, , drop = FALSE]
  for (v in cols) dfplot[[v]] <- log2(dfplot[[v]])
  list(data = dfplot, ref = 0)
}

.axis_lab <- function(name, scale) if (isTRUE(scale)) paste0("log2(", name, ")") else name
