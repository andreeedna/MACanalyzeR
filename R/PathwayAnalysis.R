PathAnalyzeR <- function(mac_obj, pathway=NULL, meta=NULL){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")

  if (is.null(mac_obj@MetaData[,meta])) {
    if (meta == "Mac") {
      stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
    } else if (meta == "Foam"){
      stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
    } else {
      stop("meta must be Sample, Cluster, Mac or Foam")
    }
  }

  if (is.null(pathway)) {
    if (mac_obj@organism=="mm") pathway <- MACanalyzeR::mouse_selected
    else pathway <- MACanalyzeR::human_selected
  } else if (!is.list(pathway)) {
    stop("Pathway must be a list! Use data('mouse_selected') or data('human_selected') to have an example")
  }

  if (is.null(meta)) meta <- mac_obj@ident

  check.cols <- function(vect){
    for (x in vect) {
      if (x < 0.001) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  path_names <- names(pathway)
  mtx <- mac_obj@CountData@NormalizedCount
  cond_tot <- mac_obj@MetaData[,meta]
  cond <- levels(factor(cond_tot))

  ll <- list()

  for (store in c("Total","SingleCell")) {
    if (store == "SingleCell") {
      path_expression <- matrix(NA,nrow=length(colnames(mtx)),ncol=length(path_names),dimnames = list(colnames(mtx), path_names))
    } else {
      path_expression <- matrix(NA,nrow=length(path_names),ncol=length(cond)+1,dimnames = list(path_names,c(cond, "pval")))
    }

    #progress bar
    message(store)
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

      if(store == "SingleCell"){
        mtx_pathway <- mtx_pathway[keep,]
        mean_expr <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
        ratio_expr <- mtx_pathway / colMeans(mean_expr)
        gene_weight <- apply(mtx_pathway, 1, var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
        path_expression[,p] <-  mean_exp_pathway
      } else {
        mtx_pathway <- mtx_pathway[keep,]
        sample_mean <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
        ratio_expr <- t(sample_mean) / colMeans(sample_mean)
        gene_weight <- apply(mtx_pathway, 1, var)
        mean_exp_pathway <- apply(ratio_expr, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
        path_expression[p, ] <-  c(mean_exp_pathway[cond], wilcox.test(ratio_expr[,1],mu=1)$p.value)
        }
    }
    close(pb)

    ll[[store]] <- as.data.frame(path_expression)
  }
  mac_obj@PathAnalyzeR[[meta]] <- ll
  return(mac_obj)
}

#####

PathHeat <- function(mac_obj, plot.by=mac_obj@ident, col=c("blue","white","red"), pval=0.05, txt.size=10){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (is.null(mac_obj@PathAnalyzeR[[plot.by]][["Total"]])) stop("Before PathHeatPlot() run the PathAnalyzeR(mac_obj, meta='",plot.by,"') function")

  dfplot <- mac_obj@PathAnalyzeR[[plot.by]][["Total"]]
  dfplot <- dfplot[which(dfplot$pval<pval),colnames(dfplot)!="pval"]

  sort_row <- c()
  for(i in colnames(dfplot)){
    select_row <- which(apply(dfplot, 1, max) == dfplot[,i])
    tmp <- rownames(dfplot)[select_row][order(dfplot[select_row,i],decreasing = T)]
    sort_row <- c(sort_row,tmp)
  }
  dfplot[is.na(dfplot)] <- 1
  b <- max((max(dfplot)-1), (1-min(dfplot)))
  mybreaks <- c(
    seq((1-b), 1-(b/3), length.out=33),
    seq(1-(b/3)+0.01, 1+(b/3), length.out=34),
    seq(1+(b/3)+0.01, 1+b,length.out=33)
  )
  color <- colorRampPalette(col)(100)
  pheatmap(dfplot[sort_row,],cluster_cols = F,cluster_rows = F,color=color,breaks = mybreaks, fontsize = txt.size, border_color = F)
}


PathDisplay <- function(mac_obj, meta=mac_obj@ident){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (is.null(mac_obj@PathAnalyzeR[[meta]][["SingleCell"]])) stop("Before PathHeatPlot() run PathAnalyzeR(mac_obj, meta='",meta,"')")

  message("Pathway Used for ", meta, " Pathway Analysis:")
  i=1
  for (p in colnames(mac_obj@PathAnalyzeR[[meta]][["SingleCell"]])) {
    message(i,'-', p)
    i=i+1
  }
}


PathCart <- function(mac_obj, pathway, plot.by=mac_obj@ident, split.by=NULL, fill=T, a=0.1, col=NULL, txt.size=15, ncol=3){
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (is.null(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])) stop("Before PathCartPlot() run PathAnalyzeR(mac_obj, meta='",plot.by,"')")

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

  all_pathway <- colnames(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])
  if(length(all_pathway)<max(pathway)) stop("Out of index, number of pathway is: ", length(all_pathway))
  pathway <- colnames(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])[pathway]

  dfplot <- cbind(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]][,pathway], mac_obj@MetaData)

  if(length(pathway) == 1) {
    plot <- ggplot(dfplot, aes(x=dfplot[,1], color=dfplot[,plot.by])) +
      labs(x=pathway[1], y="", color='', fill='') +
      geom_vline(xintercept = 1, colour = 'grey30') +
      theme(text = element_text(size = txt.size),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(),
            strip.background = element_rect(color = "black"))

    if (fill) {
      plot <- plot + geom_density(aes(fill=dfplot[,plot.by]), alpha=a)
    } else {
      plot <- plot + geom_density()
    }

    if (!is.null(col)) {
      if (length(col) == length(unique(dfplot[,plot.by]))) {
        plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
      } else {
        stop("Number of color must correspond number of ", plot.by, " : ",length(unique(dfplot[,plot.by])))
      }
    }

    if (!is.null(split.by)) {
      plot <- plot +
        facet_wrap(split.by, ncol = ncol)
    }

  } else if(length(pathway) == 2) {
    plot <- ggplot(dfplot, aes(x=dfplot[,1], y=dfplot[,2])) +
      labs(x=pathway[1], y=pathway[2], color=plot.by, fill=plot.by) +
      geom_vline(xintercept = 1, colour = 'grey30') +
      geom_hline(yintercept = 1, colour = 'grey30') +
      theme(text = element_text(size = txt.size),
            panel.background = element_blank(),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.line = element_line(),
            strip.background = element_rect(color = "black"))

    if (fill) {
      plot <- plot + stat_density_2d(aes(color = dfplot[,plot.by], fill = dfplot[,plot.by]), geom = "polygon", alpha = a, position = "identity")
    } else {
      plot <- plot + geom_density_2d(aes(color=dfplot[,plot.by]))
    }

    if (!is.null(col)) {
      if (length(col) == length(unique(dfplot[,plot.by]))) {
        plot <- plot + scale_fill_manual(values = col) + scale_color_manual(values = col)
      } else {
        stop("Number of color must correspond number of ", plot.by, " : ",length(unique(dfplot[,plot.by])))
      }
    }

    if (!is.null(split.by)) {
      plot <- plot +
        facet_wrap(split.by, ncol = ncol)
    }
  } else {
    stop("Please insert one or maximus two pathway")
  }
  return(plot)
}


PathPlot <- function(mac_obj, pathway, txt.size=14, shade="turbo", split.by=NULL, pt.size=2, ncol=2, max.cutoff=NA) {
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (is.null(mac_obj@PathAnalyzeR[[mac_obj@ident]][["SingleCell"]])) stop("Before PathHeatPlot() run PathAnalyzeR(mac_obj)")

  all_pathway <- colnames(mac_obj@PathAnalyzeR[[mac_obj@ident]][["SingleCell"]])
  if(length(all_pathway)<pathway) stop("Out of index, number of pathway is: ", length(all_pathway))

  pathway <- colnames(mac_obj@PathAnalyzeR[[mac_obj@ident]][["SingleCell"]])[pathway]

  dfplot <- cbind(mac_obj@MetaData, mac_obj@Reduction)
  dfplot$PATH <- mac_obj@PathAnalyzeR[[mac_obj@ident]][["SingleCell"]][,pathway]
  dfplot <- dfplot[order(dfplot$PATH, decreasing = F),]

  lab <- colnames(dfplot)

  if(!is.na(max.cutoff)) {
    if (max.cutoff<2) warning("It is not recommended to use a cutoff less than two")
    dfplot[which(dfplot$PATH>max.cutoff),"PATH"] <- max.cutoff
  }

  plot <- ggplot(dfplot, aes(x=dfplot[,ncol(dfplot)-2], y=dfplot[,ncol(dfplot)-1], color=PATH))+
    labs(x = lab[ncol(dfplot)-2], y=lab[ncol(dfplot)-1], color="") +
    ggtitle(pathway) +
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

###############################################

PathViolin <- function(mac_obj, pathway, plot.by=mac_obj@ident, txt.size = 15, ncol=1, col=NULL, intercept = NULL, stat=F) {
  if (class(mac_obj)!="MACanalyzeR") stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  if (is.null(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])) stop("Before PathHeatPlot() run PathAnalyzeR(mac_obj, meta='",plot.by,"')")

  all_pathway <- colnames(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])
  if(length(all_pathway)<max(pathway)) stop("Out of index, number of pathway is: ", length(all_pathway))

  pathway <- colnames(mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])[pathway]
  dfplot <- cbind(mac_obj@MetaData, mac_obj@PathAnalyzeR[[plot.by]][["SingleCell"]])
  s <- ncol(dfplot)
  plot <- list()

  for (i in pathway) {
    p1 <- eval(substitute(
      ggplot(dfplot, aes(x=dfplot[,plot.by], y=dfplot[,i], fill = dfplot[,plot.by], color = dfplot[,plot.by])) +
        geom_violin() +
        geom_boxplot(width=0.3, color="grey30", alpha=0.2) +
        labs(x="", y=i, color=plot.by, fill=plot.by) +
        theme(text = element_text(size = txt.size),
              panel.background = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(),
              strip.background = element_rect(color = "black"),
              legend.position="none")
      ,list(i = i)))

    if (!is.null(intercept)) {
      p1 <- p1 + geom_hline(yintercept = 1, colour = 'grey30')
    }

    if (!is.null(col)) {
      if (length(col) == length(unique(dfplot[,plot.by]))) {
        p1 <- p1 + scale_fill_manual(values = col) + scale_color_manual(values = col)
      } else {
        stop("Number of color must correspond number of ", plot.by, " : ", length(unique(dfplot[,plot.by])))
      }
    }

    if (stat) {
      permutazioni <- combn(unique(dfplot[,plot.by]), 2)
      for (lung in 1:ncol(permutazioni)) {
        p1 <- p1 + geom_signif(comparisons = list(as.vector(permutazioni[,lung])),
                               map_signif_level = F, color="black", margin_top = 0.05+lung*0.1)
      }
    }

    plot[[i]] <- p1  # add each plot into plot list
  }
  grid.arrange(grobs=plot, ncol=ncol)
}




