macSpectrum <- function(mac_obj, mode="macspectrum", plot.by="Sample", split.by=NULL, style="density", txt.size=17, col=NULL, ncol=3) {
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (style!="point" && style!="density") {
    stop("The style flag must be 'point' (scatter plot) or 'density' (2D density plot)")
  }

  if (is.null(mac_obj@MetaData[,plot.by])) {
    if (plot.by == "Mac") {
      stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
    } else if (plot.by == "Foam"){
      stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
    } else {
      stop("plot.by must be Sample, Cluster, Mac or Foam")
    }
  }

  if (!is.null(split.by)) {
    if (is.null(mac_obj@MetaData[,split.by])) {
      if (split.by == "Mac") {
        stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
      } else if (split.by == "Foam"){
        stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
      } else {
        stop("split.by must be Sample, Cluster, Mac or Foam")
      }
    }
  }

  macspec <- function(mac_mtx, org) {
    if (org!="mm" && org!="hs") {
      stop("The org flag must be 'mm' (Mus musculus) or 'hs' (Homo sapiens)")
    }

    MFDI <- function(mac_mtx, org) {
      con_mean <- MACanalyzeR::con_mean_new
      foam_mean <- MACanalyzeR::foam_mean_new

      id <- "symbol"
      if (org!="mm" || id!="symbol") {
        dict <- MACanalyzeR::human_mouse_dict

        names(con_mean) <- dict[names(con_mean),paste(org,id, sep="_")]
        con_mean <- con_mean[!is.na(names(con_mean))]
        names(foam_mean) <- dict[names(foam_mean),paste(org,id, sep="_")]
        foam_mean <- foam_mean[!is.na(names(foam_mean))]
      }

      MFI_genes <- intersect(names(foam_mean), rownames(mac_mtx))
      con_mean <- con_mean[MFI_genes]
      foam_mean <- foam_mean[MFI_genes]
      MFI_mtx <- mac_mtx[MFI_genes,]

      #sigma of mac cells
      mac_sigma <- 1 : ncol(MFI_mtx)
      total_gene_number <- nrow(MFI_mtx)
      for(i in 1:ncol(MFI_mtx)) {
        mac_sigma[i]<-(sum(MFI_mtx[,i]^2)/total_gene_number)^0.5
      }

      #sigma of con mean:
      total_gene_number <- length(con_mean)
      con_sigma <- (sum(con_mean^2)/total_gene_number)^0.5

      #sigma of foam mean:
      total_gene_number <- length(foam_mean)
      foam_sigma <- (sum(foam_mean^2)/total_gene_number)^0.5

      #correlation of mac - con mean:
      total_gene_number <- length(foam_mean)
      con_pearson<-1:ncol(MFI_mtx)
      for (j in 1:ncol(MFI_mtx)){
        con_pearson[j]<-sum((MFI_mtx[,j]/mac_sigma[j])*(con_mean/con_sigma))/total_gene_number
      }

      #correlation of mac - foam mean:
      total_gene_number<-length(foam_mean)
      foam_pearson<-1:ncol(MFI_mtx)
      for (j in 1:ncol(MFI_mtx)){
        foam_pearson[j]<-sum((MFI_mtx[,j]/mac_sigma[j])*(foam_mean/foam_sigma))/total_gene_number
      }

      a<-1
      b<-1
      c<- 0
      x0<-foam_pearson
      y0<-con_pearson
      d_sqr<-(a*x0+b*y0+c)^2/(a^2+b^2)
      x_start<--1
      y_start<-(-a)*x_start+(-c)
      x_end<-1
      y_end<-(-a)*x_end+(-c)

      l<-((x0-x_start)^2+(y0-y_start)^2-d_sqr)^0.5
      l_max<-((x_end-x_start)^2+(y_end-y_start)^2-d_sqr)^0.5
      MFDI<-(l-0)/(l_max-0)*100-50

      return(MFDI)
    }

    M1_mean <- MACanalyzeR::M1_mean_new
    M2_mean <- MACanalyzeR::M2_mean_new
    M0_mean <- MACanalyzeR::M0_mean_new

    id <- "symbol"
    if (org!="mm" || id!="symbol") {
      dict <- MACanalyzeR::human_mouse_dict

      names(M1_mean) <- dict[names(M1_mean),paste(org,id, sep="_")]
      M1_mean <- M1_mean[!is.na(names(M1_mean))]

      names(M2_mean) <- dict[names(M2_mean),paste(org,id, sep="_")]
      M2_mean <- M2_mean[!is.na(names(M2_mean))]

      names(M0_mean) <- dict[names(M0_mean),paste(org,id, sep="_")]
      M0_mean <- M0_mean[!is.na(names(M0_mean))]
    }

    mac_mtx <- mac_mtx-rowMeans(mac_mtx)

    ##  MPI matrix
    MPI_genes <- intersect(names(M1_mean),rownames(mac_mtx))
    M1_mean <- M1_mean[MPI_genes]
    M2_mean <- M2_mean[MPI_genes]
    MPI_mtx <- mac_mtx[MPI_genes,]

    ## ADMI matrix
    AMDI_genes <- intersect(names(M0_mean),rownames(mac_mtx))
    M0_mean <- M0_mean[AMDI_genes]
    AMDI_mtx <- mac_mtx[AMDI_genes,]

    #sigma of MPI
    mac_sigma <- 1:ncol(MPI_mtx)
    total_gene_number<-nrow(MPI_mtx)
    for(i in 1:ncol(MPI_mtx)) {
      #options(digits=9)
      mac_sigma[i]<-(sum(MPI_mtx[,i]^2)/total_gene_number)^0.5
    }

    #sigma of ADMI:
    mac_sigma_m0<-1:ncol(AMDI_mtx)
    total_gene_number<-nrow(AMDI_mtx)
    for(i in 1:ncol(AMDI_mtx)) {
      #options(digits=9)
      mac_sigma_m0[i]<-(sum(AMDI_mtx[,i]^2)/total_gene_number)^0.5
    }

    #sigma of M0 mean:
    total_gene_number <- length(M0_mean)
    M0_sigma <- (sum(M0_mean^2)/total_gene_number)^0.5

    #sigma of M1 mean:
    total_gene_number <- length(M1_mean)
    M1_sigma <- (sum(M1_mean^2)/total_gene_number)^0.5

    #sigma of M2 mean:
    total_gene_number <- length(M2_mean)
    M2_sigma <- (sum(M2_mean^2)/total_gene_number)^0.5


    #correlation of mac - M0 mean:
    total_gene_number <- nrow(AMDI_mtx)
    M0_pearson <- 1:ncol(AMDI_mtx)
    for (j in 1:ncol(AMDI_mtx)){
      M0_pearson[j]<-sum((AMDI_mtx[,j]/mac_sigma_m0[j])*(M0_mean/M0_sigma))/total_gene_number
    }

    #correlation of Mac - M1 mean:
    total_gene_number <- length(M2_mean)
    M1_pearson<-1:ncol(MPI_mtx)
    for (j in 1:ncol(MPI_mtx)){
      M1_pearson[j]<-sum((MPI_mtx[,j]/mac_sigma[j])*(M1_mean/M1_sigma))/total_gene_number
    }

    #correlation of mac - M2 mean:
    total_gene_number <- length(M2_mean)
    M2_pearson<-1:ncol(MPI_mtx)
    for (j in 1:ncol(MPI_mtx)){
      M2_pearson[j]<-sum((MPI_mtx[,j]/mac_sigma[j])*(M2_mean/M2_sigma))/total_gene_number
    }

    a <-0.991414467
    b <-1
    c <- -0.0185412856

    d_sqr<-(a*M1_pearson+b*M2_pearson+c)^2/(a^2+b^2)
    x_start<--1
    y_start<-(-a)*x_start+(-c)
    x_end<-1
    y_end<-(-a)*x_end+(-c)

    l<-((M1_pearson-x_start)^2+(M2_pearson-y_start)^2-d_sqr)^0.5
    l_max<-((x_end-x_start)^2+(y_end-y_start)^2-d_sqr)^0.5
    MPI<-(l-0)/(l_max-0)*100-50

    AMDI<- -M0_pearson*50

    ms <- data.frame("MPI" = MPI,
                     "AMDI" = AMDI,
                     "MDFI" = MFDI(mac_mtx, org),
                     row.names=colnames(mac_mtx),
                     stringsAsFactors=F)

    return(ms)
  }

  mac_obj@MetaData <- cbind(mac_obj@MetaData, macspec(mac_obj@CountData@NormalizedCount, org = mac_obj@organism))
  dfplot <- mac_obj@MetaData

  if (mode == "macspectrum") {
    plot <- ggplot(dfplot, aes(x=MPI, y=AMDI)) +
      labs(fill=plot.by, color=plot.by) +
      geom_vline(xintercept = 1, colour = 'grey30') +
      geom_hline(yintercept = 0.5, colour = 'grey30') +
      theme(text = element_text(size = txt.size),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(),
            strip.background = element_rect(color = "black"))
  } else if (mode == "atherospectrum") {
    plot <- ggplot(dfplot, aes(x=MPI, y=MDFI)) +
      labs(fill=plot.by, color=plot.by) +
      geom_vline(xintercept = 1, colour = 'grey30') +
      geom_hline(yintercept = 0.5, colour = 'grey30') +
      theme(text = element_text(size = txt.size),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(),
            strip.background = element_rect(color = "black"))
  } else {
    stop("The mode flag must be 'MACanalyzeR' or 'atherospectrum'")
  }

  if (style == 'point') {
    plot <- plot +
      geom_point(aes(color=dfplot[,plot.by]), size=2)
  } else {
    plot <- plot +
      stat_density_2d(aes(color = dfplot[,plot.by], fill = dfplot[,plot.by]),
                      geom = "polygon", alpha = 0.2, position = "identity")
  }

  if (!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  if (!is.null(col) && length(col)== length(unique(dfplot[,plot.by]))) {
    plot <- plot +
      scale_color_manual(values = col) +
      scale_fill_manual(values = col)
  } else if (!is.null(col)) {
    warning("Number of color must correspond number of plot.by classes")
  }

  return(plot)
}
