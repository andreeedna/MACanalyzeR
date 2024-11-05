FoamSpotteR <- function(mac_obj){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (mac_obj@organism=="mm") FoamModel <- MACanalyzeR::mouse_LAMprey
  else FoamModel <- MACanalyzeR::human_LAMprey

  FoamGene <- rownames(FoamModel$importance)
  TOTgene <- intersect(FoamGene, mac_obj@CountData@Feature)
  Foamtx <- t(mac_obj@CountData@NormalizedCount[TOTgene,])

  if (!identical(TOTgene, FoamGene)) {
    stop("All the genes for the prediction are not present in the NormalizedCountData")
  }


  FoamProb <- as.data.frame(predict(FoamModel, Foamtx, type = "prob", decision.value = TRUE))
  FoamPred <- predict(FoamModel, Foamtx)
  FoamProb$Foam <- FoamPred

  mac_obj@MetaData <- cbind(mac_obj@MetaData, FoamProb)

  return(mac_obj)
}


