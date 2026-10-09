\name{CreateMacObj}
\alias{CreateMacObj}

\title{Create a MACanalyzeR object from a Seurat object}

\description{
Builds the \code{MACanalyzeR} object used by every function of the package, starting from a preprocessed and clustered \code{Seurat} object (Seurat v4 and v5 are both supported). The function extracts the normalized expression matrix, the cell embeddings of the chosen dimensionality reduction and the cell metadata.
}

\usage{
CreateMacObj(seurat_obj, id = "seurat_clusters", reduction = "umap", org = "mm",
             assay = NULL, sparse = TRUE)
}

\arguments{
  \item{seurat_obj}{A \code{Seurat} object that has already been normalized, reduced and clustered. For best results it should contain only monocytes/macrophages.}
  \item{id}{Name of the column of \code{seurat_obj@meta.data} holding the main cell grouping (clusters or cell annotations). It becomes the default grouping (\code{mac_obj@ident}) for the \code{plot.by} and \code{meta} arguments of the other functions. Default: \code{"seurat_clusters"}.}
  \item{reduction}{Name of the dimensionality reduction stored in \code{seurat_obj@reductions} whose embeddings are used for plotting (e.g. \code{"umap"}, \code{"tsne"}). Only the first two dimensions are kept. Default: \code{"umap"}.}
  \item{org}{Organism of the dataset: \code{"mm"} (\emph{Mus musculus}) or \code{"hs"} (\emph{Homo sapiens}). It selects the organism-specific models and gene sets used by the analysis modules. Default: \code{"mm"}.}
  \item{assay}{Name of the assay whose normalized data are used (e.g. \code{"RNA"} or \code{"SCT"}). Default: \code{NULL}, the default assay of the Seurat object.}
  \item{sparse}{Logical. If \code{TRUE}, a sparse expression matrix is kept sparse, which greatly reduces memory usage; each module only densifies the genes it needs. Use \code{FALSE} to store a dense matrix as in previous versions. Default: \code{TRUE}.}
}

\details{
The normalized counts are taken from the \code{data} layer of the selected assay. With Seurat v5, if the \code{data} layer is split (e.g. after \code{split()} or integration), join the layers first with \code{SeuratObject::JoinLayers()}. Subsetting the Seurat object to the myeloid compartment beforehand is strongly recommended.

All the columns of \code{seurat_obj@meta.data} are copied, with their original names, into \code{mac_obj@MetaData}. Any of them can therefore be passed to the \code{plot.by}, \code{split.by} and \code{meta} arguments of the other functions. The analysis modules add further columns: \code{Foam}, \code{fMAC+} and \code{fMAC-} (\code{\link{FoamSpotteR}}) and \code{Mac} (\code{\link{MacPolarizeR}}).
}

\value{
An object of class \code{MACanalyzeR} with the following slots:
  \item{CountData}{Normalized expression matrix (\code{@NormalizedCount}, genes x cells, sparse or dense) and gene names (\code{@Feature}).}
  \item{Reduction}{Data frame with the first two dimensions of the cell embeddings.}
  \item{MetaData}{Data frame with the cell metadata.}
  \item{MacPolarizeR, PathAnalyzeR, MitoScanneR}{Empty lists, filled by the corresponding modules.}
  \item{ident}{Name of the default grouping column (the \code{id} argument).}
  \item{organism}{The \code{org} argument.}
  \item{version}{Package version used to create the object.}
}

\seealso{
\code{\link{MacPlot}}, \code{\link{FoamSpotteR}}, \code{\link{MacPolarizeR}}, \code{\link{PathAnalyzeR}}, \code{\link{MitoScanneR}}
}

\examples{
\dontrun{
seurat_obj <- readRDS("macrophages.rds")
head(seurat_obj@meta.data)

mac_obj <- CreateMacObj(seurat_obj,
                        id = "Cells",
                        reduction = "umap",
                        org = "mm")
}
}
