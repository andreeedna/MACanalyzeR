MacPolarizeR <- function(mac_obj){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (mac_obj@organism=="mm") inf <- MACanalyzeR::mouse_PolGenes
  else inf <- MACanalyzeR::human_PolGenes

  mac_mtx <- mac_obj@CountData@NormalizedCount
  m1 <- intersect(inf$M1, rownames(mac_mtx))
  m2 <- intersect(inf$M2, rownames(mac_mtx))
  genes <- c(m1,m2)

  m_mtx <- t(as.data.frame(mac_mtx[genes,]))
  mac_km <- kmeans(m_mtx, 3, iter.max = 5000)

  # annotation
  for (store in c("Total","SingleCell")) {
    if (store == "SingleCell") {
      path_expression <- matrix(NA, nrow=length(colnames(mac_mtx)), ncol=2, dimnames=list(colnames(mac_mtx), c("M1", "M2")))
    } else {
      path_expression <- matrix(NA, nrow=2, ncol=3, dimnames=list(c("M1", "M2"),c("1","2","3")))
    }
    for (p in c("M1", "M2")) {
      genes <- inf[[p]]
      genes <- intersect(genes, rownames(mac_mtx))
      mtx_pathway <- mac_mtx[genes,]
      mtx_pathway <- mtx_pathway[rowSums(mtx_pathway)>0,]
      sample_mean <- apply(mtx_pathway, 1, function(x)by(x, mac_km$cluster, mean))

      if(store == "SingleCell"){
        mean_expr <- apply(mtx_pathway, 1, function(x)by(x, mac_km$cluster, mean))
        ratio_expr <- mtx_pathway / colMeans(mean_expr)
        gene_weight <- apply(mtx_pathway, 1, var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
        path_expression[,p] <-  mean_exp_pathway
      } else {
        sample_mean <- apply(mtx_pathway, 1, function(x)by(x, mac_km$cluster, mean))
        ratio_expr <- t(sample_mean) / colMeans(sample_mean)
        gene_weight <- apply(mtx_pathway, 1, var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
        path_expression[p, ] <- mean_exp_pathway
      }
    }
    mac_obj@MacPolarizeR[[store]] <- path_expression
  }

  mi <- which.max(mac_obj@MacPolarizeR$Total["M1",])
  li <- which.max(mac_obj@MacPolarizeR$Total["M2",])
  if (mi == li) {
    colnames(mac_obj@MacPolarizeR$Total)[mi] <- "Inflammatory"
    li <- which.max(mac_obj@MacPolarizeR$Total["M2",-mi])
    colnames(mac_obj@MacPolarizeR$Total)[li] <- "Healing"
    colnames(mac_obj@MacPolarizeR$Total)[-c(mi,li)] <- "Transitional"
  } else {
    colnames(mac_obj@MacPolarizeR$Total)[mi] <- "Inflammatory"
    colnames(mac_obj@MacPolarizeR$Total)[li] <- "Healing"
    colnames(mac_obj@MacPolarizeR$Total)[-c(mi,li)] <- "Transitional"
  }


  mac_obj@MetaData$Mac <- factor(mac_km$cluster, levels = 1:3, labels = colnames(mac_obj@MacPolarizeR$Total))
  mac_obj@MetaData$Mac <- factor(mac_obj@MetaData$Mac, levels = c("Inflammatory", "Transitional", "Healing"))

  return(mac_obj)
}

