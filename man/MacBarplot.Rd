\name{MacBarplot}
\alias{MacBarplot}

\title{
Macrophage polarization Barplot
}

\description{
The MacBarplot() function creates a barplot for visualizing polarization classes proportion.
}

\usage{
MacBarplot(mac_obj,
           plot.by = "Sample",
           txt.size = 14,
           col = c("#C5283D", "#E9724C", "#FFC857")
           )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mac_obj}{
MACanalyzer object
}
  \item{plot.by}{
Name of one metadata to group the MacBarplot (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{txt.size}{
Text size of MacBarplot [default: 14]
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted (3).
}
}

\examples{
##  MacPolarizeR prediction
mac_obj <- MacPolarizeR(mac_obj)

##  BarPlot
MacBarplot(mac_obj)
}
