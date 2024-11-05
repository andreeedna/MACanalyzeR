\name{FoamPlot}
\alias{FoamPlot}

\title{
FoamDEX Plot
}

\description{
The FoamPlot() function creates a 2D scatter plot for visualizing FoamDEX score on single-cell.
Each point in the plot represents a cell coloured by the FoamDEX score.
}

\usage{
FoamPlot(mac_obj,
         txt.size = 14,
         shade = "viridis",
         split.by = NULL,
         pt.size = 2,
         ncol = 2)
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{txt.size}{
Text size of FoamPlot [default: 14]
}
  \item{shade}{
Shades of FoamDEX plot. Eight option are available: magma, inferno, plasma, viridis, cividis, rocket, mako and turbo [default: viridis]
}
  \item{split.by}{
Name of one metadata to split the FoamPlot by (Sample, Cluster, Mac, Foam) [default: NULL]
}
  \item{pt.size}{
Point size of FoamPlot [default: 2]
}
  \item{ncol}{
Number of columns for the FoamPlot split
}
}

\details{
The FoamDEX is a value between 0 and 1 that indicated the probability of that cell is a foam-like macrophage
}

\examples{
##  FoamSpotteR prediction
mac_obj <- FoamSpotteR(mac_obj)

##  base FoamPlot
FoamPlot(mac_obj)

##  splitter FoamPlot
FoamPlot(mac_obj, split.by="Sample")
}
