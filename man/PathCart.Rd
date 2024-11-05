\name{PathCart}
\alias{PathCart}

\title{
PathAnalyzeR Cartplot
}
\description{
PathAnalyzeR implement a novel and unconventional visualization approach called Cartplot, which utilizes 2D density plots. The `PathCart()` function allows to simultaneously observe the expression distribution of two pathways at the single-cell level. This unique visualization provides valuable insights into the coordinated upregulation of pathways within the context of the analyzed metadata.
}

\usage{
PathCart(mac_obj,
         pathway,
         plot.by = "Sample",
         fill = F,
         a = 0.1,
         col = NULL,
         txt.size = 15
         )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{pathway}{
A vector with pathways' number to visualize [MAX 2]
}
  \item{plot.by}{
Name of one metadata to group the cells by (Sample, Cluster, Mac, Foam) [default: Sample]
}
  \item{fill}{
Logical. Shows the plot filled
}
  \item{a}{
Alpha value of the filled FoamLine [default: 0.5]
}
  \item{col}{
A vector of colors. The number of colors must match the number of classes to be plotted.
}
  \item{txt.size}{
Text size of FoamLine [default: 14]
}
}

\note{
Is available also single pathway density plot, specifying one number in the pathway flag
}

\examples{
# PathAnalyzeR pathway score
mac_obj <- PathAnalyzeR(mac_obj)

# Pathway Cartplot
PathCart(mac_obj, pathway=c(1,2))
}
