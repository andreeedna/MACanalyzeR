\name{MitoBalance}
\alias{MitoBalance}

\title{
Mitochondrial Balance gene expression ratio
}
\description{
The `MitoBalance()` function generates a violin plot to visualize potential shifts in the expression levels of mitochondrial and mitonuclear genes. This plot depicts the ratio between the expression of these two gene types.
}

\usage{
MitoBalance(mac_obj,
            plot.by = "Sample",
            txt.size = 15,
            ncol = 1,
            col = NULL,
            intercept = F,
            title = "MitoNuclear Balance"
            )
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
Number of columns for the MitoBalance plot
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
MitoBalance(mac_obj)
}
