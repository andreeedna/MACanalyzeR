FoamPlot <- function(mac_obj, txt.size=14, shade="viridis", split.by=NULL, pt.size=2, ncol=3){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (is.null(mac_obj@MetaData$`fMAC+`)){
    stop("FoamSpotteR prediction is absent. Please, before run FoamSpotteR() function")
  }

  if (!is.null(split.by)) {
    if (!split.by %in% colnames(mac_obj@MetaData)) {
      if (split.by == "Mac") {
        stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
      } else {
        stop("split.by must be in colnames(mac_obj@MetaData)")
      }
    }
  }

  dfplot <- cbind(mac_obj@MetaData, mac_obj@Reduction)
  dfplot <- dfplot[order(dfplot$`fMAC+`, decreasing = F),]
  lab <- colnames(dfplot)

  plot <- ggplot(dfplot, aes(x=dfplot[,ncol(dfplot)-1], y=dfplot[,ncol(dfplot)], color=`fMAC+`))+
    labs(x = lab[ncol(dfplot)-1], y=lab[ncol(dfplot)], color="FoamDEX") +
    geom_point(size=pt.size) +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black")) +
    scale_color_viridis_c(option = shade)

  if(!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  return(plot)
}


###################################

FoamLine <- function(mac_obj, plot.by=mac_obj@ident, split.by=NULL, txt.size=14, fill=T, a=0.5, col=NULL, ncol=3){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (is.null(mac_obj@MetaData$`fMAC+`)){
    stop("FoamSpotteR prediction is absent. Please, before run FoamSpotteR() function")
  }

  if (!plot.by %in% colnames(mac_obj@MetaData)) {
    if (plot.by == "Mac") {
      stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
    } else {
      stop("plot.by must be in colnames(mac_obj@MetaData)")
    }
  }

  if(!is.null(split.by)){
    if (is.null(mac_obj@MetaData[,split.by])) {
      if (split.by == "Mac") {
        stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
      } else if (split.by == "Foam"){
        stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
      } else {
        stop("split.by must be in colnames(mac_obj@MetaData")
      }
    }
  }

  dfplot <- mac_obj@MetaData
  lab <- colnames(dfplot)

  plot <- ggplot(dfplot, aes(x=`fMAC+`, color=dfplot[,plot.by]))+
    labs(x = "FoamDEX", y='', fill=plot.by, color=plot.by) +
    xlim(0,1) +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"))

  if (fill) {
    plot <- plot + geom_density(aes(y = after_stat(scaled), fill=dfplot[,plot.by]), alpha=a)
  } else {
    plot <- plot + geom_density(aes(y = after_stat(scaled)))
  }

  if (!is.null(col)) {
    if (length(col) == length(unique(dfplot[,plot.by]))) {
      plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
    } else {
      stop("Number of color must correspond number of ", plot.by, " : ",length(unique(dfplot[,1])))
    }
  }

  if(!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  return(plot)
}



