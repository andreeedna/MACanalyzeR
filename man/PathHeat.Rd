\name{PathHeat}
\alias{PathHeat}

\title{
PathAnalyzeR Heatmap
}

\description{
The PathHeat() function plots the result of PathAnalyzeR metadata analysis as an heatmap that allows a general view of the enriched pathways between the specified group.
}

\usage{
PathHeat(mac_obj,
         plot.by = "Sample",
         col = c("blue", "white", "red"),
         pval = 0.05,
         txt.size=10
         )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{plot.by}{
Name of one metadata to group the PathHeat (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{col}{
a vector of color, used as a shade for the PathHeat.
}
  \item{pval}{
pvalue threshold [default: 0.05]
}
  \item{txt.size}{
Text size of PathHeat [default: 10]
}
}

\note{
Before plotting PathHeat for a specific metadata, you must launch PathAnalyzeR() function with the appropriate meta flag corresponding to that metadata
}

\examples{
mac_obj <- PathAnalyzer(mac_obj, meta="Sample")
PathHeat(mac_obj, meta="Sample")

mac_obj <- PathAnalyzer(mac_obj, meta="Cluster")
PathHeat(mac_obj, meta="Cluster")
}
