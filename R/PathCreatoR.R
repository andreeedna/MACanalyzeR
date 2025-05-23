PathCreatoR <- function(organism, subcollection, as.df = FALSE) {
  # Check if organism is valid
  if (!(organism %in% c("hs", "mm"))) {
    stop("Organism not found. Please select one of: hs, mm")
  }

  # Check if subcollection is valid
  if (organism == "hs") {
    gsea_data <- gsea_human
  } else {
    gsea_data <- gsea_mouse
  }
  if (!(subcollection %in% names(gsea_data))) {
    stop("Subcollection not found. Please select one of: ", paste(names(gsea_data), collapse = ", "))
  }

  # Extract pathway data
  pathway_data <- gsea_data[[subcollection]]

  # Convert to data frame if requested
  if (as.df) {
    pathway_df <- do.call(rbind, lapply(names(pathway_data), function(term) {
      data.frame(TERM = term, Gene = pathway_data[[term]], stringsAsFactors = FALSE)
    }))
    return(pathway_df)
  }

  # Return pathway data as list
  return(pathway_data)
}