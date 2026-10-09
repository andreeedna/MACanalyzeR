setClass("Macmatrices",
         slots = list(NormalizedCount = "ANY",
                      Feature = "vector"))

setClass("MACanalyzeR",
         slots = list(CountData = "Macmatrices",
                      Reduction = "data.frame",
                      MetaData = "data.frame",
                      MacPolarizeR = "list",
                      PathAnalyzeR = "list",
                      MitoScanneR = "list",
                      ident = "vector",
                      organism = "vector",
                      version = "vector"))

setMethod("show", "MACanalyzeR", function(object) {
  m <- object@CountData@NormalizedCount
  cat("An object of class MACanalyzeR (v", as.character(object@version), ")\n", sep = "")
  cat(nrow(m), "features across", ncol(m), "cells -",
      if (identical(object@organism, "hs")) "Homo sapiens" else "Mus musculus", "\n")
  cat("Default grouping:", object@ident, "\n")
  done <- c(FoamSpotteR = "Foam" %in% colnames(object@MetaData),
            MacPolarizeR = "Mac" %in% colnames(object@MetaData))
  cat("FoamSpotteR:", if (done["FoamSpotteR"]) "done" else "not run", "| MacPolarizeR:",
      if (done["MacPolarizeR"]) "done" else "not run", "\n")
  if (length(object@PathAnalyzeR) > 0) cat("PathAnalyzeR:", paste(names(object@PathAnalyzeR), collapse = ", "), "\n")
  if (length(object@MitoScanneR) > 0) cat("MitoScanneR:", paste(names(object@MitoScanneR), collapse = ", "), "\n")
  invisible(object)
})


## Normalized expression from a Seurat object (v4 or v5).
.seurat_data <- function(seurat_obj, assay) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found. Available assays: ",
         paste(names(seurat_obj@assays), collapse = ", "), call. = FALSE)
  }
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    a <- seurat_obj[[assay]]
    if (methods::is(a, "Assay5")) {
      lyr <- SeuratObject::Layers(a, search = "data")
      if (length(lyr) == 0) {
        stop("Assay '", assay, "' has no 'data' layer. Run Seurat::NormalizeData() first.", call. = FALSE)
      }
      if (length(lyr) > 1) {
        stop("The 'data' layer of assay '", assay, "' is split into ", length(lyr), " layers (",
             paste(lyr, collapse = ", "), "). Run seurat_obj <- SeuratObject::JoinLayers(seurat_obj) first.",
             call. = FALSE)
      }
    }
    SeuratObject::LayerData(seurat_obj, assay = assay, layer = "data")
  } else {
    SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = "data")
  }
}


CreateMacObj <- function(seurat_obj, id = "seurat_clusters", reduction = "umap", org = "mm",
                         assay = NULL, sparse = TRUE) {
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("The input object is not a Seurat object", call. = FALSE)
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required to read Seurat objects", call. = FALSE)
  }
  if (!org %in% c("mm", "hs")) {
    stop("The org argument must be 'mm' (Mus musculus) or 'hs' (Homo sapiens)", call. = FALSE)
  }
  if (!id %in% colnames(seurat_obj@meta.data)) {
    stop("The column '", id, "' is not present in seurat_obj@meta.data", call. = FALSE)
  }
  if (is.null(seurat_obj@reductions[[reduction]])) {
    stop("The reduction '", reduction, "' is not present in the Seurat object. Available reductions: ",
         paste(names(seurat_obj@reductions), collapse = ", "), call. = FALSE)
  }
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat_obj)

  pkg_version <- as.character(utils::packageVersion("MACanalyzeR"))
  message("MACanalyzeR version: ", pkg_version)
  message("Organism: ", if (org == "mm") "Mus musculus" else "Homo sapiens")
  message("Seurat object version: ", as.character(seurat_obj@version), " | assay: ", assay)

  mat <- .seurat_data(seurat_obj, assay)
  if (!sparse || !methods::is(mat, "sparseMatrix")) mat <- as.matrix(mat)

  cells <- colnames(mat)
  emb <- seurat_obj@reductions[[reduction]]@cell.embeddings
  if (ncol(emb) < 2) stop("The reduction '", reduction, "' has less than two dimensions", call. = FALSE)
  if (!all(cells %in% rownames(emb))) stop("Some cells are missing from the reduction '", reduction, "'", call. = FALSE)
  emb <- as.data.frame(emb[cells, 1:2, drop = FALSE])

  meta <- as.data.frame(seurat_obj@meta.data)[cells, , drop = FALSE]

  if (ncol(mat) > 20000 && !methods::is(mat, "sparseMatrix")) {
    message("Note: storing a dense matrix for ", ncol(mat), " cells may require a lot of memory")
  }

  new("MACanalyzeR",
      CountData = new("Macmatrices", NormalizedCount = mat, Feature = rownames(mat)),
      Reduction = emb,
      MetaData = meta,
      MacPolarizeR = list(),
      PathAnalyzeR = list(),
      MitoScanneR = list(),
      ident = id,
      organism = org,
      version = pkg_version)
}



MacPlot <- function(mac_obj, plot.by = mac_obj@ident, pt.size = 2, txt.size = 14,
                    split.by = NULL, ncol = 3, col = NULL) {
  .check_obj(mac_obj)
  .check_col(mac_obj, plot.by, "plot.by")
  .check_col(mac_obj, split.by, "split.by")

  emb <- .embedding(mac_obj)
  dfplot <- cbind(mac_obj@MetaData, emb$data)

  plot <- ggplot(dfplot, aes(x = .data[[".dim1"]], y = .data[[".dim2"]], color = .data[[plot.by]])) +
    geom_point(size = pt.size) +
    labs(x = emb$labels[1], y = emb$labels[2], color = "") +
    .mac_theme(txt.size)

  if (!is.null(split.by)) plot <- plot + facet_wrap(split.by, ncol = ncol)

  if (plot.by == "Mac" && is.null(col)) {
    col <- c("Inflammatory" = "#C5283D", "Transitional" = "#E9724C", "Healing" = "#FFC857")
  }
  if (!is.null(col)) {
    if (length(col) == length(unique(stats::na.omit(dfplot[[plot.by]])))) {
      plot <- plot + scale_color_manual(values = col)
    } else {
      warning("The number of colours must match the number of groups of ", plot.by, " (",
              length(unique(stats::na.omit(dfplot[[plot.by]]))), "): using the default palette", call. = FALSE)
    }
  }

  plot
}
