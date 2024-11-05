\name{MitoScanneR}
\alias{MitoScanneR}

\title{
Mitochondrial metabolism Scanner
}

\description{
The MitoScanneR module is dedicated to the analysis of oxidative and mitochondrial profile. This module allows the identification of cells with high oxidative and mitochondrial acrivity. It's based on PathAnalyzeR pathway score on MitoCarta genes.
}

\usage{
mac_obj <- MitoScanneR(mac_obj,
                       meta = "Sample"
                       )
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{meta}{
metadata selected to perform MitoScanneR
}
}

\note{
This version of MACanalyzeR uses genes from MitoCarta 3.0
}

\examples{
mac_obj <- MitoScanneR(mac_obj)
}
