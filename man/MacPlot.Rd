\name{MacPlot}
\alias{MacPlot}

\title{Plot cells on the dimensionality reduction}

\description{
Draws a 2D scatter plot of the cells on the embedding stored in the \code{MACanalyzeR} object (e.g. UMAP or t-SNE), coloured by a metadata column and optionally split into panels.
}

\usage{
MacPlot(mac_obj, plot.by = mac_obj@ident, pt.size = 2, txt.size = 14,
        split.by = NULL, ncol = 3, col = NULL)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{plot.by}{Name of the metadata column used to colour the cells. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}. Use \code{"Foam"} or \code{"Mac"} to show the results of \code{\link{FoamSpotteR}} or \code{\link{MacPolarizeR}}.}
  \item{pt.size}{Point size. Default: \code{2}.}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. If its length does not match the number of groups a warning is raised and the default palette is used. When \code{plot.by = "Mac"} and no colours are given, the MacPolarizeR palette is used.}
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{FoamPlot}}, \code{\link{PathPlot}}}

\examples{
\dontrun{
MacPlot(mac_obj)

## colour by sample and split by cluster
MacPlot(mac_obj, plot.by = "Sample", split.by = "Cluster")

## show the MacPolarizeR classes
mac_obj <- MacPolarizeR(mac_obj)
MacPlot(mac_obj, plot.by = "Mac")
}
}
