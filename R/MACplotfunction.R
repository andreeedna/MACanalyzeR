MacBarplot <- function(mac_obj, plot.by=mac_obj@ident, txt.size=14, col=c("#C5283D", "#E9724C", "#FFC857")){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (is.null(mac_obj@MetaData[,"Mac"])){
    stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
  }

  if (is.null(mac_obj@MetaData[,plot.by])) {
    if (plot.by == "Foam"){
      stop("FoamSpotteR prediction is absent. Please, before run FoamSpotteR() function")
    } else {
      stop("plot.by must be in colnames(mac_obj@MetaData)")
    }
  }

  count <- table(mac_obj@MetaData[,plot.by], mac_obj@MetaData$Mac)
  count <- apply(count, 2, function(x) x*100/rowSums(count))

  dfplot <- merge(factor(rownames(count)), factor(c("Inflammatory", "Transitional", "Healing"), levels = c("Inflammatory", "Transitional", "Healing")))
  #dfplot <- merge(rownames(count), colnames(count))
  dfplot$value <- c(count)

  plot <- ggplot(dfplot, aes(x=x, y=value, fill=y)) +
    geom_bar(stat = "identity") +
    labs(x="", y="", fill="") +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black", size=txt.size),
          axis.text.y = element_text(color = "black", size=txt.size),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"))

  if(length(col) == 3){
    plot <- plot + scale_fill_manual(values=col)
  } else {
    message("The vector must contain 3 colors!")
    plot <- plot + scale_fill_manual(values=c("Inflammatory" = "#C5283D", "Transitional" = "#E9724C", "Healing" = "#FFC857"))
  }

  return(plot)
}



MacRadar <- function(mac_obj, plot.by=mac_obj@ident, split.by=NULL, sort=F, txt.size=14, col=c("#C5283D", "#E9724C", "#FFC857"), ncol=3){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (is.null(mac_obj@MetaData[,"Mac"])){
    stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
  }

  if (is.null(mac_obj@MetaData[,plot.by])) {
    if (plot.by == "Foam"){
      stop("FoamSpotteR prediction is absent. Please, before run FoamSpotteR() function")
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

  dfplot <- mac_obj@MetaData
  if (sort) {
    dfplot[,plot.by] <- factor(dfplot[,plot.by], levels = rev(names(sort(table(dfplot[,plot.by])))), ordered = T)
  } else {
    dfplot[,plot.by] <- factor(dfplot[,plot.by])
  }

  plot <- ggplot(dfplot, aes(x=dfplot[,plot.by], fill=Mac)) +
    geom_bar() +
    coord_polar(start = -(pi/length(levels(dfplot[,plot.by])))) +
    labs(fill="") +
    ylim(c(-(max(table(dfplot[,plot.by])))/5, max(table(dfplot[,plot.by])))) +
    theme_void() +
    theme(text = element_text(size = txt.size),
          axis.text.x = element_text(colour = "black"))

  if (!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  if(length(col) == 3){
    plot <- plot + scale_fill_manual(values=col)
  } else {
    message("The vector must contain 3 colors!")
    plot <- plot + scale_fill_manual(values=c("Inflammatory" = "#C5283D", "Transitional" = "#E9724C", "Healing" = "#FFC857"))
  }

  return(plot)
}


PolCart <- function(mac_obj, plot.by=mac_obj@ident, split.by=NULL, style="density", fill=T, a=0.1, col=NULL, txt.size=15, ncol=3){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")

  if (is.null(mac_obj@MetaData[,"Mac"])){
    stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
  }

  if (is.null(mac_obj@MetaData[,plot.by])) {
    if (plot.by == "Foam"){
      stop("FoamSpotteR prediction is absent. Please, before run FoamSpotteR() function")
    } else {
      stop("plot.by must be Sample, Cluster, Mac or Foam")
    }
  }

  if(!is.null(split.by)){
    if (!split.by %in% colnames(mac_obj@MetaData)) {
      if (split.by == "Mac") {
        stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
      } else if (split.by == "Foam"){
        stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
      } else {
        stop("split.by must be in colnames(mac_obj@MetaData")
      }
    }
  }

  dfplot <- cbind(mac_obj@MetaData, mac_obj@MacPolarizeR$SingleCell)

  plot <- ggplot(dfplot, aes(x=dfplot[,ncol(dfplot)-1], y=dfplot[,ncol(dfplot)], color=dfplot[,plot.by], fill=dfplot[,plot.by])) +
    labs(x="M1", y="M2", color=plot.by, fill=plot.by) +
    geom_vline(xintercept = 1, colour = 'grey30') +
    geom_hline(yintercept = 1, colour = 'grey30') +
    theme(text = element_text(size = 17),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"))

  if (!is.null(col)) {
    plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
  }

  if (style == 'point') {
    plot <- plot +
      geom_point(size=2)
  } else {
    plot <- plot +
      stat_density_2d(geom = "polygon", alpha = 0.2, position = "identity")
  }

  if (!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  return(plot)
}




