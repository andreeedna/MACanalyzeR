\name{MacPolarizeR}
\alias{MacPolarizeR}

\title{
Macrophage Polarization Analysis
}

\description{
The MacPolarizeR module classify monocytes and macrophages based on their polarization phenotype. This module operates by clustering cells according to the expression of genes associated with M1 and M2 polarization and cluster cells into 3 different class: Inflammatory (M1-like), Healing (M2-like) and Transitional (M0-like).

In addition to to classification, MacPolarizeR calculates an M1 and M2 polarization score for each cell based on the expression of the gene markers. This score system facilitates a more comprehensive and nuanced characterization of macrophages polarization state
}

\usage{
mac_obj <- MacPolarizeR(mac_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mac_obj}{
MACanalyzer object
}
}

\note{
The output of the cell polarization clustering are stored within the metadata of the MACanalyzeR object accessible at mac_obj@MetaData$Mac and the polarization score is accessible at mac_obj@MacPolarizeR
}

\examples{
##  MacPolarizeR classification
mac_obj <- MacPolarizeR(mac_obj)
}
