setClass("Macmatrices",
         slots = list(NormalizedCount='matrix',
                      Feature="vector"
                      )
         )

setClass("MACanalyzeR",
         slots = list(CountData="Macmatrices",
                      Reduction="data.frame",
                      MetaData="data.frame",
                      MacPolarizeR='list',
                      PathAnalyzeR="list",
                      MitoScanneR='list',
                      ident='vector',
                      organism='vector',
                      version='vector'
         )
)

# setMethod("show","MACanalyzeR",function(object){
#   p <- m <- c()
#   for (w in c("Sample", "Cluster", "Mac", "Foam")) {
#     p <- append(p, w)
#     p <- append(p, if(!length(slot(object@PathAnalyzeR, w))==0) "- Calculated |"  else "- Not Calculated |")
#
#     m <- append(m, w)
#     m <- append(m, if(!length(slot(object@MitoScanneR, w))==0) "- Calculated |"  else "- Not Calculated |")
#   }
#
#   cat("An object of class MACanalyzeR\n",
#       nrow(object@CountData@NormalizedCount), " features across ", ncol(object@CountData@NormalizedCount), " cells", " | ",
#       length(unique(object@MetaData$Sample)), " Sample(s) - ", length(unique(object@MetaData$Cluster)), " Cluster(s)", "\n",
#       "\n",
#       # "MacPolarizeR:", "\t",
#       # if (is.null(object@MetaData[["Mac"]])) "Not Calculated" else "Calculated", "\n",
#       # "FoamSpotteR:", "\t",
#       # if (is.null(object@MetaData[["Foam"]])) "Not Calculated" else "Calculated", "\n",
#       # "PathAnalyzeR:", "\t",
#       # paste(p, collapse=" "), "\n",
#       # "MitoScanneR:", "\t",
#       # paste(m, collapse=" "), "\n",
#       sep='')
# }
# )
#
# setMethod("print","MACanalyzeR", show())

CreateMacObj <- function(seurat_obj, id="seurat_clusters", reduction = 'umap', org="mm"){
  message("MACanalyzeR version: 1.0.1")

  if (org == "mm") {
    message("Organism: Mus musculus")
  } else if (org == "hs"){
    message("Organism: Homo sapiens")
  } else {
    stop("The organism flag should be mm (Mus musculus) or hs (Homo sapiens)")
  }

  if(class(seurat_obj) == "Seurat"){
    #qua deve esserci un controllo se è oggetto 4 o oggetto 5
    if (as.integer(substr(seurat_obj@version, start = 1, stop = 1)) < 5) {
      #oggetto 4
      message("Seurat version: ", seurat_obj@version)
      mat_norm <- as.matrix(seurat_obj@assays$RNA@data)
    } else if (as.integer(substr(seurat_obj@version, start = 1, stop = 1)) >= 5) {
      #oggetto 5
      message("Seurat version: ", seurat_obj@version)
      mat_norm <- as.matrix(seurat_obj@assays$RNA$data)
    }

    if(id %in% colnames(seurat_obj@meta.data)){
        if (!is.null(seurat_obj@reductions[[reduction]])) {
          mat <- new("Macmatrices",
                     NormalizedCount=mat_norm,
                     Feature=rownames(mat_norm))

          mac_obj <- new("MACanalyzeR",
                         CountData=mat,
                         Reduction=as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings),
                         MetaData=as.data.frame(seurat_obj@meta.data),
                         MacPolarizeR=list(),
                         PathAnalyzeR=list(),
                         MitoScanneR=list(),
                         ident=id,
                         organism=org,
                         version='1.0.0'
          )

          return(mac_obj)
        } else {
          stop("This reduction is not present in the seurat object. Please insert a valid reduction [try 'umap' or 'tsne']")
        }
      } else {
      stop("The column ", id, " is not present in seurat_obj@metadata")
    }
  } else {
    stop("The input object is not a SeuratObject. Please insert a SeuratObject")
  }
}



MacPlot <- function(mac_obj, plot.by=mac_obj@ident, pt.size=2 ,txt.size=14, split.by=NULL, ncol=3, col=NULL){
  if (class(mac_obj)!="MACanalyzeR"){
    stop("The input object is not a MACanalyzeR Object. Please insert a MACanalyzeR Object.
  Create it with CreateMacObj() function.")
  }

  if (!plot.by %in% colnames(mac_obj@MetaData)) {
    if (plot.by == "Mac") {
      stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
    } else if (plot.by == "Foam"){
      stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
    } else {
      stop("plot.by must be in colnames(mac_obj@MetaData)")
    }
  }

  if (!is.null(split.by)) {
    if (!split.by %in% colnames(mac_obj@MetaData)) {
      if (split.by == "Mac") {
        stop("MacPolarizeR prediction is absent. Please, before run MacPolarizeR() function")
      } else if (split.by == "Foam"){
        stop("FoamFinder prediction is absent. Please, before run FoamSpotteR() function")
      } else {
        stop("split.by must be in colnames(mac_obj@MetaData)")
      }
    }
  }

  dfplot <- cbind(mac_obj@MetaData, mac_obj@Reduction)
  labs <- colnames(dfplot)

  plot <- ggplot(dfplot, aes(x=dfplot[,ncol(dfplot)-1], y=dfplot[,ncol(dfplot)], color=dfplot[,plot.by])) +
    geom_point(size=pt.size) +
    labs(x=labs[ncol(dfplot)-1], y=labs[ncol(dfplot)], color="") +
    theme(text = element_text(size = txt.size),
          panel.background = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.line = element_line(),
          strip.background = element_rect(color = "black"))

  if (!is.null(split.by)) {
    plot <- plot +
      facet_wrap(split.by, ncol = ncol)
  }

  if (plot.by == "Mac") col = c("Inflammatory" = "#C5283D", "Transitional" = "#E9724C", "Healing" = "#FFC857")
  if (!is.null(col)) {
    if (length(col) == length(unique(mac_obj@MetaData[,plot.by]))) {
      plot <- plot + scale_color_manual(values = col)
    } else {
      warning("Number of color must correspond number of ", plot.by, " : ",length(unique(mac_obj@MetaData[,plot.by])))
    }
  }

  return(plot)
}






