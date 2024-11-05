\name{PathCreatoR}
\alias{PathCreatoR}

\title{
Pathway list Creator
}
\description{
PathCreatoR function download gene sets from msigdb in a MACanalyzeR-compatible structure.
}
\usage{
p <- PathCreatoR(organism, collections = NULL, subcollections = NULL)
}

\arguments{
  \item{organism}{
Reference organism to download the gene sets ["mm" or "hs"]
}
  \item{collections}{
MSIGDB collection of interest
}
  \item{subcollections}{
MSIGDB subcollection of interest
}
}
\details{

}

\examples{
p <- PathCreatoR("mm", collections="c2")
}
