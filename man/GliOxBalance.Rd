\name{GliOxBalance}
\alias{GliOxBalance}

\title{
Glicolitic Oxidative Balance
}
\description{
The `GliOxBalance()` function determine whether there is a shift in the expression of genes related to glycolytic or oxidative phosphorylation.
}
\usage{
GliOxBalance(mac_obj,
            plot.by = "Cluster",
            txt.size = 15,
            ncol = 1,
            col = NULL,
            intercept = F,
            title = "GlicoOxphos Balance")
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{plot.by}{
Name of one metadata to group the cells by (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{txt.size}{
Text size of FoamLine [default: 14]
}
  \item{ncol}{
Number of columns for the GliOxBalance plot
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted.
}
  \item{intercept}{
Logical. Plot the intercept y=1
}
  \item{title}{
Plot title
}
}

\examples{
GliOxBalance(mac_obj)

}
