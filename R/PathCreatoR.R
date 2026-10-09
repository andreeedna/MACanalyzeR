#' Fetch MSigDB Collections or Convert Custom Gene Lists
#'
#' @param organism "hs" for human, "mm" for mouse
#' @param subcollection MSigDB collection name
#' @param as.df Logical; return a data frame for clusterProfiler enricher() or GSEA()
#' @param gene.list Optional; a custom list of gene sets to convert to DF
#' @return A list of gene sets usable for PathAnalyzeR or a data frame for clusterProfiler enricher() or GSEA()
#' @export
PathCreatoR <- function(
  organism = "hs",
  subcollection = NULL,
  as.df = FALSE,
  gene.list = NULL
) {
  # Check if gene.list is provided
  if (!is.null(gene.list)) {
    # Warning if both gene.list and collection are provided
    if (!is.null(subcollection)) {
      warning(
        "Both subcollection and gene.list provided. Ignoring subcollection and using custom gene.list"
      )
    }

    # Stop if gene.list is a list of character vectors
    if (!is.list(gene.list)) {
      stop("gene.list must be a list of character vectors")
    }

    # If provided, convert gene.list to a tidy data frame
    if (as.df) {
      df <- stack(gene.list)
      colnames(df) <- c("Gene", "TERM")
      df$Gene <- as.character(df$Gene)
      df$TERM <- as.character(df$TERM)
      return(df[, c("TERM", "Gene")])
    }

    # Return the gene.list as is
    return(gene.list)
  }

  # Interface with MSigDB
  if (!(organism %in% c("hs", "mm"))) {
    stop("Organism not found. Please select one of: 'hs', 'mm'")
  }

  # Dynamically assign data based on organism
  # These datasets are assumed to be available in the package data
  if (organism == "hs") {
    gsea_data <- gsea_human
  } else {
    gsea_data <- gsea_mouse
  }

  # Check if subcollection is valid
  if (is.null(subcollection) || !(subcollection %in% names(gsea_data))) {
    stop(
      "Subcollection not found. Please select one of: ",
      paste(names(gsea_data), collapse = ", ")
    )
  }

  # Fetch MSigDB collection data
  pathway_data <- gsea_data[[subcollection]]

  # If provided, convert pathway_data to a tidy data frame
  if (as.df) {
    pathway_df <- stack(pathway_data)
    colnames(pathway_df) <- c("Gene", "TERM")
    pathway_df$Gene <- as.character(pathway_df$Gene)
    pathway_df$TERM <- as.character(pathway_df$TERM)
    return(pathway_df[, c("TERM", "Gene")])
  }

  # Return the pathway_data as is
  return(pathway_data)
}
