\name{CreateMacObj}
\alias{CreateMacObj}

\title{
Create MACanalyzeR object
}

\description{
Create the MACanalyzer object from a Seurat object (compatible with Seurat v4 and v5). The function extract information from Seurat object metadata corresponding to the Sample, the experimental design (treated/not treated), and to the Cluster, the macrophage subclustering.
}

\usage{
mac_obj <- CreateMacObj(seurat_obj,
                        sample,
                        cluster = "seurat_clusters",
                        reduction = "umap"
                        )
}

\arguments{
  \item{seurat_obj}{
Seurat object already processed and clusterized
}
  \item{sample}{
Column name in seurat_obj@meta.data corresponding to the experimental design
}
  \item{cluster}{
Column name in seurat_obj@meta.data corresponding to cluster information or cell annotation (usually used for macrophage subclusters) [default: seurat_clusters]
}
  \item{reduction}{
Name of redution utilized for cell visualization [default: umap]
}
}

\note{
The MACanalyzeR object will have four main metadata:

-Sample: The experimental design

-Cluster: The macrophage subtypes (or any other annotation)

-Foam: The FoamSpotteR() classification

-Mac: The MacPolarizeR() classification

In the MACanalyzeR functions, the `plot.by` and `meta` flags take one of these four arguments as input (the name used in the Seurat object is not preserved!)
}

\examples{
seurat_obj <- readRDS("vignette/mac.rds")
head(seurat_obj@meta.data)

mac_obj <- CreateMacObj(seurat_obj,
                        "Sample",
                        "Cells",
                        reduction="tsne")
}

