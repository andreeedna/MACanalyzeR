\name{MacRadar}
\alias{MacRadar}

\title{
Macrophage polarization Radarplot
}
\description{
The MacRadar() function creates a radarplot for visualizing polarization classes proportion.
}

\usage{
MacRadar(mac_obj,
         plot.by = "Sample",
         txt.size = 14,
         col = c("#C5283D", "#E9724C", "#FFC857")
         )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{plot.by}{
Name of one metadata to group the MacRadar (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{txt.size}{
Text size of MacRadar [default: 14]
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted (3).
}
}


\examples{
##  MacPolarizeR prediction
mac_obj <- MacPolarizeR(mac_obj)

##  BarPlot
MacRadar(mac_obj)
}
