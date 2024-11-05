\name{PathViolin}
\alias{PathViolin}

\title{
Pathway Violinplot
}

\description{
 The `PathViolin()` function allows you to visualize the distribution of expression for individual pathways and enables statistical comparisons between these expressions across different metadata.
 }

\usage{
PathViolin(mac_obj,
           pathway,
           plot.by = "Sample",
           txt.size = 15,
           ncol = 1,
           col = NULL,
           intercept = NULL,
           stat = F
           )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{pathway}{
A vector with pathways' number to visualize
}
  \item{plot.by}{
Name of one metadata to group the cells by (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{txt.size}{
Text size of FoamLine [default: 14]
}
  \item{ncol}{
Number of columns for the PathViolin
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted.
}
  \item{intercept}{
Logical. Plot the intercept y=1
}
  \item{stat}{
Logical. Show the pvalues
}
}

\examples{
# PathAnalyzeR pathway score
mac_obj <- PathAnalyzeR(mac_obj)

# Pathway Violinplot
PathViolin(mac_obj, pathway=c(1,2))
}
