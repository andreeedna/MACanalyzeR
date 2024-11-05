\name{FoamLine}
\alias{FoamLine}

\title{
FoamDEX density Line
}

\description{
The FoamLine() function creates a density plot of FoamDEX for individual samples or clusters.
This graphical representation allows users to easily compare how FoamDEX levels are spread across different conditions.
}

\usage{
FoamLine(mac_obj,
         plot.by = "Sample",
         txt.size = 14,
         fill = F,
         a = 0.5,
         col = NULL
         )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{plot.by}{
Name of one metadata to group (color) cells by (Sample, Cluster, Mac, Foam) [default: Cluster]
}
  \item{txt.size}{
Text size of FoamLine [default: 14]
}
  \item{fill}{
Logical.
}
  \item{a}{
Alpha value of the filled FoamLine [default: 0.5]
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted.
}
}

\examples{
##  FoamSpotteR prediction
mac_obj <- FoamSpotteR(mac_obj)

##  base FoamLine
FoamLine(mac_obj)

##  splitter FoamPlot
FoamPlot(mac_obj,
         plot.by="Cluster",
         fill=T)
}
