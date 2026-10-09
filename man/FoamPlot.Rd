\name{FoamPlot}
\alias{FoamPlot}

\title{Plot the FoamDEX score on the dimensionality reduction}

\description{
Draws a 2D scatter plot of the cells on the embedding, coloured by their FoamDEX score (the probability of being a foam-like macrophage computed by \code{\link{FoamSpotteR}}). Cells with the highest scores are drawn on top.
}

\usage{
FoamPlot(mac_obj, txt.size = 14, shade = "viridis", split.by = NULL,
         pt.size = 2, ncol = 3)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{FoamSpotteR}} has been run.}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{shade}{Viridis colour scale used for the FoamDEX. One of \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"viridis"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"} or \code{"turbo"}. Default: \code{"viridis"}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{pt.size}{Point size. Default: \code{2}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
}

\details{
The FoamDEX ranges from 0 to 1 and represents the probability that a cell is a foam-like macrophage. The colour scale always spans the full 0-1 range, so plots of different datasets are directly comparable.
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{FoamSpotteR}}, \code{\link{FoamLine}}}

\examples{
\dontrun{
mac_obj <- FoamSpotteR(mac_obj)

FoamPlot(mac_obj)
FoamPlot(mac_obj, split.by = "Sample", shade = "magma", pt.size = 1.5)
}
}
