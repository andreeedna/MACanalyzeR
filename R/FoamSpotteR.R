FoamSpotteR <- function(mac_obj, max.missing = 0.1) {
  .check_obj(mac_obj)

  FoamModel <- .org_data(mac_obj, MACanalyzeR::mouse_LAMprey, MACanalyzeR::human_LAMprey)
  FoamGene <- rownames(FoamModel$importance)
  missing <- setdiff(FoamGene, rownames(mac_obj@CountData@NormalizedCount))

  if (length(missing) / length(FoamGene) > max.missing) {
    stop(length(missing), " of the ", length(FoamGene), " genes of the FoamSpotteR model are missing from the data (",
         paste(missing, collapse = ", "), "). Use the full gene set, not only variable features.", call. = FALSE)
  }

  Foamtx <- matrix(0, nrow = ncol(mac_obj@CountData@NormalizedCount), ncol = length(FoamGene),
                   dimnames = list(colnames(mac_obj@CountData@NormalizedCount), FoamGene))
  present <- setdiff(FoamGene, missing)
  Foamtx[, present] <- t(.expr(mac_obj, present))[, present, drop = FALSE]
  if (length(missing) > 0) {
    warning(length(missing), " model gene(s) missing from the data were set to 0: ",
            paste(missing, collapse = ", "), ". Predictions may underestimate foam-like cells.", call. = FALSE)
  }
  Foamtx <- as.data.frame(Foamtx, check.names = FALSE)

  FoamProb <- predict(FoamModel, Foamtx, type = "prob")
  FoamPred <- predict(FoamModel, Foamtx)

  ## assign the columns (re-running FoamSpotteR overwrites them instead of duplicating them)
  mac_obj@MetaData[["fMAC-"]] <- as.numeric(FoamProb[, "fMAC-"])
  mac_obj@MetaData[["fMAC+"]] <- as.numeric(FoamProb[, "fMAC+"])
  mac_obj@MetaData[["Foam"]] <- factor(as.character(FoamPred), levels = c("fMAC-", "fMAC+"))

  tab <- table(mac_obj@MetaData$Foam)
  message("FoamSpotteR: ", tab[["fMAC+"]], " fMAC+ and ", tab[["fMAC-"]], " fMAC- cells")
  mac_obj
}
