#' Create MACanalyzeR-compatible Seurat Object
#'
#' @param object A Seurat object
#' @param group.by Column name in meta.data to use as identity (e.g. "seurat_clusters")
#' @param organism Organism: "mm" (mouse) or "hs" (human)
#' @return A Seurat object with MACanalyzeR metadata initialized
#' @export
CreateMacObj <- function(object, group.by = "seurat_clusters", organism = "mm") {
  message("MACanalyzeR version: 2.0.0")

  if (organism == "mm") {
    message("Organism: Mus musculus")
  } else if (organism == "hs") {
    message("Organism: Homo sapiens")
  } else {
    stop("Organism must be 'mm' (Mus musculus) or 'hs' (Homo sapiens)")
  }

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!(group.by %in% colnames(object@meta.data))) {
    stop("Metadata column '", group.by, "' not found")
  }

  # Initialize MACanalyzeR storage in @misc
  if (is.null(object@misc$MACanalyzeR)) {
    object@misc$MACanalyzeR <- list(
      organism = organism,
      ident = group.by,
      version = "2.0.0"
    )
  }

  return(object)
}

