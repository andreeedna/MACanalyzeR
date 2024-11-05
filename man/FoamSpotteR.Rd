\name{FoamSpotteR}
\alias{FoamSpotteR}

\title{
Foam-like macrophage Spotter
}

\description{
The FoamSpotteR module enables the identification of cell exhibiting foam-like characteristic. This function employs a randomForest binary classifier on the data, assigning a label (fMAC+ for foam-like macrophages and fMAC- for non foam-like macrophage) to each cell.
Addictionally, it computes a score, FoamDEX, for each cell, representing the probability of being a foam-like macrophage.
}

\usage{
mac_obj <- FoamSpotteR(mac_obj)
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
}

\note{
The details regarding the FoamDEX score and the prediction result are stored in the MACanalyzeR object at mac_obj@FoamSpotteR
}

\examples{
##  FoamSpotteR prediction
mac_obj <- FoamSpotteR(mac_obj)
head(mac_obj@FoamSpotteR)
}
