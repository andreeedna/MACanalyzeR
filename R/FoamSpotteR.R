#' Foam Cell Identification (FoamSpotteR)
#'
#' @param object A Seurat object
#' @param assay Assay to use (default: DefaultAssay(object))
#' @param layer Layer to use (default: "data")
#' @return A Seurat object with foam cell probabilities
#' @importFrom stats predict
#' @importFrom Seurat GetAssayData DefaultAssay
#' @export
FoamSpotteR <- function(object, assay = NULL, layer = "data") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  organism <- mac_info$organism
  if (organism == "mm") {
    FoamModel <- MACanalyzeR::mouse_LAMprey
  } else {
    FoamModel <- MACanalyzeR::human_LAMprey
  }

  FoamGene <- rownames(FoamModel$importance)
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)
  TOTgene <- intersect(FoamGene, rownames(mtx))

  if (!identical(sort(TOTgene), sort(FoamGene))) {
    stop("Not all genes required for prediction are present in the dataset.")
  }

  Foamtx <- t(as.matrix(mtx[FoamGene, ]))

  FoamProb <- as.data.frame(stats::predict(FoamModel, Foamtx, type = "prob"))
  FoamPred <- stats::predict(FoamModel, Foamtx)
  FoamProb$FoamLabel <- FoamPred

  # Add to meta.data
  colnames(FoamProb) <- paste0("Foam_", colnames(FoamProb))
  object@meta.data <- cbind(object@meta.data, FoamProb)

  return(object)
}


