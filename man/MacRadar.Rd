\name{MacRadar}
\alias{MacRadar}

\title{Radar plot of macrophage polarization classes}

\description{
Draws a circular (radar-like) barplot in which each group is a slice and the bar height is the number of cells, coloured by the polarization class assigned by \code{\link{MacPolarizeR}}. It gives an overview of both the size of each group and its polarization composition.
}

\usage{
MacRadar(mac_obj, plot.by = mac_obj@ident, split.by = NULL, sort = FALSE,
         txt.size = 14, col = c("#C5283D", "#E9724C", "#FFC857"), ncol = 3)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{MacPolarizeR}} has been run.}
  \item{plot.by}{Name of the metadata column defining the groups (one slice per group). Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{sort}{Logical. If \code{TRUE}, groups are ordered by decreasing number of cells. Default: \code{FALSE}.}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{col}{Vector of three colours for the Inflammatory, Transitional and Healing classes, in this order. If a different number of colours is given, the default palette is used.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{MacPolarizeR}}, \code{\link{MacBarplot}}}

\examples{
\dontrun{
mac_obj <- MacPolarizeR(mac_obj)

MacRadar(mac_obj, plot.by = "Cluster", sort = TRUE)
MacRadar(mac_obj, plot.by = "Cluster", split.by = "Sample")
}
}
