\name{PathPlot}
\alias{PathPlot}

\title{
Pathway MacPlot
}
\description{
Through the `PathPlot()` function, we can visualize the enrichment score for each cell independently within a MacPlot. This enables the visualization of pathway enrichment patterns across the entire cell embedding.
}
\usage{
PathPlot(mac_obj,
         pathway,
         txt.size = 14,
         shade = "turbo",
         split.by = NULL,
         pt.size = 2,
         ncol = 2,
         max.cutoff = NA)
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{pathway}{
A vector with pathways' number to visualize
}
  \item{txt.size}{
Text size of FoamLine [default: 14]
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
  \item{max.cutoff}{
Maximum value of pathway score
}
}

\examples{
# PathAnalyzeR pathway score
mac_obj <- PathAnalyzeR(mac_obj)

# Pathway MacPlot
PathPlot(mac_obj, pathway=c(1,2))
}
