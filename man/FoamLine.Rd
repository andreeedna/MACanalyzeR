\name{FoamLine}
\alias{FoamLine}

\title{FoamDEX density plot}

\description{
Draws the distribution of the FoamDEX score (computed by \code{\link{FoamSpotteR}}) as one density curve per group, making it easy to compare foam-like features across samples, clusters or any other metadata.
}

\usage{
FoamLine(mac_obj, plot.by = mac_obj@ident, split.by = NULL, txt.size = 14,
         fill = TRUE, a = 0.5, col = NULL, ncol = 3)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{FoamSpotteR}} has been run.}
  \item{plot.by}{Name of the metadata column defining the groups (one curve per group). Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{fill}{Logical. If \code{TRUE} the area under each curve is filled. Default: \code{TRUE}.}
  \item{a}{Transparency (alpha) of the filled areas, between 0 and 1. Used only when \code{fill = TRUE}. Default: \code{0.5}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
}

\details{
Densities are scaled so that each curve has a maximum of 1, which makes groups with very different numbers of cells directly comparable. The x axis spans the whole FoamDEX range (0 to 1).
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{FoamSpotteR}}, \code{\link{FoamPlot}}}

\examples{
\dontrun{
mac_obj <- FoamSpotteR(mac_obj)

FoamLine(mac_obj)
FoamLine(mac_obj, plot.by = "Sample", fill = FALSE)
FoamLine(mac_obj, plot.by = "Sample", split.by = "Cluster")
}
}
