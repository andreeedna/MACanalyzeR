MitoScanneR <- function(mac_obj, meta=mac@ident){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")

  if (is.null(mac_obj@MetaData[[meta]])) {
    if (meta == "MAC") {
      stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
    } else if (meta == "Foam"){
      stop("LAMFinder prediction is absent. Please, before run LAMFindeR() function")
    } else {
      stop("meta must be Sample, Cluster, MAC or Foam")
    }
  }

  check.cols <- function(vect){
    for (x in vect) {
      if (x < 0.001) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  if (mac_obj@organism=="mm") MitoPath <- MACanalyzeR::mouse_MitoCarta
  else MitoPath <- MACanalyzeR::human_MitoCarta

  MitoList <- list()

  mtx <- mac_obj@CountData@NormalizedCount
  cond_tot <- mac_obj@MetaData[[meta]]
  cond <- unique(cond_tot)

  for (m in names(MitoPath)) {
    pathway <- MitoPath[[m]]
    path_names <- names(pathway)
    mito_expression <- matrix(NA,nrow=length(path_names),ncol=length(cond),dimnames = list(path_names,cond))

    #progress bar
    message(m)
    pb <- txtProgressBar(min = 0, max = length(path_names)-1, style = 3, width = 60, char = "+")
    i=0

    for (p in path_names) {
      ####### progress bar ###
      setTxtProgressBar(pb, i)
      i = i + 1
      ########################

      genes <- pathway[[p]]
      genes <- intersect(genes, rownames(mtx))
      if(length(genes) < 5) next
      mtx_pathway <- mtx[genes,]
      mtx_pathway <- mtx_pathway[rowSums(mtx_pathway)>0,]
      sample_mean <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
      keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
      if(length(keep)<3) next
      mtx_pathway <- mtx_pathway[keep,]
      mean_expr <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
      ratio_expr <- t(mean_expr) / colMeans(mean_expr)
      gene_weight <- apply(mtx_pathway, 1, var)
      mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
      mito_expression[p,] <-  mean_exp_pathway
    }
    close(pb)

    MitoList[[m]] <- as.data.frame(mito_expression)
  }
  mac_obj@MitoScanneR[[meta]] <- MitoList

  return(mac_obj)
}

HeatMito <- function(mac_obj, plot.by=mac@ident, col=c("blue","white","red"), txt.size=10){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (length(mac_obj@MitoScanneR[[plot.by]])==0) stop("Before HeatMito() run the MitoScanneR(mac_obj, meta='",plot.by,"') function")

  for (mitoname in names(mac_obj@MitoScanneR[[plot.by]][8:1])) {
    dat <- mac_obj@MitoScanneR[[plot.by]][[mitoname]]
    sort_row <- c()
    for(i in colnames(dat)){
      select_row <- which(apply(dat, 1, max) == dat[,i])
      tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
      sort_row <- c(sort_row,tmp)
    }
    dat[is.na(dat)] <- 1
    b <- max((max(dat)-1), (1-min(dat)))
    mybreaks <- c(
      seq((1-b), 1-(b/3), length.out=33),
      seq(1-(b/3)+0.01, 1+(b/3), length.out=34),
      seq(1+(b/3)+0.01, 1+b,length.out=33)
    )
    color <- colorRampPalette(col)(100)
    pheatmap(dat[sort_row,],cluster_cols = F,cluster_rows = F,color=color,breaks = mybreaks, main=gsub("_"," ",mitoname),
             fontsize = txt.size, border_color=F)
  }

}

MitoBalance <- function(mac_obj, plot.by = "Sample", txt.size = 15, ncol=1, col=NULL, intercept = F, title="MitoNuclear Balance") {
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")

  check.cols <- function(vect){
    for (x in vect) {
      if (x < 0.001) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  mtx <- mac_obj@CountData@NormalizedCount
  cond_tot <- as.vector(mac_obj@MetaData[[plot.by]])
  cond <- unique(cond_tot)
  path_names <- c("OXPHOS mitochondrial subunits", "OXPHOS nuclear subunits")

  path_expression <- matrix(NA,nrow=length(colnames(mtx)),ncol=length(path_names),dimnames = list(colnames(mtx), path_names))

  for (p in path_names) {
    ##  genes <- MACanalyzeR::MitoScan$OXPHOS[[p]]
    if (mac_obj@organism=="mm") genes <- MACanalyzeR::mouse_MitoCarta$OXPHOS[[p]]
    else genes <- MACanalyzeR::human_MitoCarta$OXPHOS[[p]]

    genes <- intersect(genes, rownames(mtx))
    if(length(genes) < 5) next
    mtx_pathway <- mtx[genes,]
    mtx_pathway <- mtx_pathway[rowSums(mtx_pathway)>0,]
    sample_mean <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if(length(keep)<3) next

    mtx_pathway <- mtx_pathway[keep,]
    mean_expr <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
    ratio_expr <- mtx_pathway / colMeans(mean_expr)
    gene_weight <- apply(mtx_pathway, 1, var)
    mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
    path_expression[,p] <-  mean_exp_pathway
  }

  dat <- as.data.frame(path_expression)
  dat[,"MitochondrialBalance"] <- log10(dat[,1]/dat[,2])
  dat[,plot.by] <- mac_obj@MetaData[[plot.by]]

  threshold <- !is.infinite(dat$MitochondrialBalance)
  message("Removed ", sum(!threshold), " cell containing non-finite value")
  dat <- dat[threshold,]

  plot <- list()

  p1 <- ggplot(dat, aes(x=dat[,4], y=dat[,3], fill = dat[,4], color = dat[,4])) +
    geom_violin() +
    geom_boxplot(width=0.3, color="grey30", alpha=0.2) +
    labs(x="", y="Mitochondrial/Nuclear GeneRatio") +
    ggtitle(title) +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          legend.key=element_rect(fill="white"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"),
          legend.position="none")

  if (intercept) {
    p1 <- p1 + geom_hline(yintercept = 0, colour = 'grey30')
  }

  if (!is.null(col)) {
    if (length(col) == length(unique(dat[,plot.by]))) {
      p1 <- p1 + scale_fill_manual(values = col) + scale_color_manual(values = col)
    } else {
      stop("Number of color must correspond number of ", plot.by, " : ", length(unique(dat[,4])))
    }
  }
  return(p1)
}

GliOxBalance <- function(mac_obj, plot.by = "Cluster", txt.size = 15, ncol=1, col=NULL, intercept = F, title="GlicoOxphos Balance") {
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")

  check.cols <- function(vect){
    for (x in vect) {
      if (x < 0.001) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  mtx <- mac_obj@CountData@NormalizedCount
  cond_tot <- as.vector(mac_obj@MetaData[[plot.by]])
  cond <- unique(cond_tot)
  path_names <- c("GLYCOLYSIS_GLUCONEOGENESIS", "OXIDATIVE_PHOSPHORYLATION")

  path_expression <- matrix(NA,nrow=length(colnames(mtx)),ncol=length(path_names),dimnames = list(colnames(mtx), path_names))

  for (p in path_names) {
    if (mac_obj@organism=="mm") genes <- MACanalyzeR::mouse_selected[[p]]
    else genes <- MACanalyzeR::human_selected[[p]]

    genes <- intersect(genes, rownames(mtx))
    if(length(genes) < 5) next
    mtx_pathway <- mtx[genes,]
    mtx_pathway <- mtx_pathway[rowSums(mtx_pathway)>0,]
    sample_mean <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
    keep <- colnames(sample_mean)[apply(sample_mean, 2, check.cols)]
    if(length(keep)<3) next

    mtx_pathway <- mtx_pathway[keep,]
    mean_expr <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
    ratio_expr <- mtx_pathway / colMeans(mean_expr)
    gene_weight <- apply(mtx_pathway, 1, var)
    mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
    path_expression[,p] <-  mean_exp_pathway
  }

  dat <- as.data.frame(path_expression)
  dat[,"MitochondrialBalance"] <- log10(dat[,1]/dat[,2])
  dat[,plot.by] <- mac_obj@MetaData[[plot.by]]

  threshold <- !is.infinite(dat$MitochondrialBalance)
  message("Removed ", sum(!threshold), " cell containing non-finite value")
  dat <- dat[threshold,]

  plot <- list()

  p1 <- ggplot(dat, aes(x=dat[,4], y=dat[,3], fill = dat[,4], color = dat[,4])) +
    geom_violin() +
    geom_boxplot(width=0.3, color="grey30", alpha=0.2) +
    labs(x="", y="Glicolysis/OXPHOS GeneRatio") +
    ggtitle(title) +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          legend.key=element_rect(fill="white"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"),
          legend.position="none")

  if (intercept) {
    p1 <- p1 + geom_hline(yintercept = 0, colour = 'grey30')
  }

  if (!is.null(col)) {
    if (length(col) == length(unique(dat[,plot.by]))) {
      p1 <- p1 + scale_fill_manual(values = col) + scale_color_manual(values = col)
    } else {
      stop("Number of color must correspond number of ", plot.by, " : ", length(unique(dat[,4])))
    }
  }
  return(p1)
}

