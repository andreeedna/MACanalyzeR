PathCreatoR <- function(organism, subcollection){
  if (organism == "hs"){
    if (!(subcollection %in% names(gsea_human))) {
      stop("Subcollection not found. Please select one of: ", paste(names(gsea_human), collapse = ", "))
    }
    return (gsea_human[[subcollection]])
  }

  else if (organism == "mm"){
    if (!(subcollection %in% names(gsea_mouse))) {
      stop("Subcollection not found. Please select one of: ", paste(names(gsea_mouse), collapse = "\n"))
    }
    return (gsea_mouse[[subcollection]])
  }

  else{
    stop("Organism not found. Please select one of: hs, mm")
  }
}

