\name{MacPlot}
\alias{MacPlot}

\title{
MACanalyzeR Plot
}

\description{
The MacPlot() function creates a 2D scatter plot for visualizing single-cell.
Each point in the plot represents a cell, positioned based on its embedding coordinates obtained from a dimensionality reduction technique.
}

\usage{
MacPlot(mac_obj,
        plot.by = "Cluster",
        pt.size = 2,
        txt.size = 14,
        split.by = NULL,
        ncol = 3,
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
  \item{pt.size}{
Point size of MacPlot [default: 2]
}
  \item{txt.size}{
Text size of MacPlot [default: 14]
}
  \item{split.by}{
Name of one metadata to split the MacPlot by (Sample, Cluster, Mac, Foam) [default: NULL]
}
  \item{ncol}{
Number of columns for the MacPlot split
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted.
}
}

\examples{
## base MacPlot
MacPlot(mac_obj)

# splitted MacPlot
MacPlot(mac_obj, plot.by="Sample", split.by="Cluster")
}
