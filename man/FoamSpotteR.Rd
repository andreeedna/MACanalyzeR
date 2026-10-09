\name{FoamSpotteR}
\alias{FoamSpotteR}

\title{Identify foam-like macrophages}

\description{
The FoamSpotteR module identifies cells with a foam-like phenotype. A random forest binary classifier, trained on foamy and non-foamy macrophages from atherosclerotic plaques, labels each cell as \code{fMAC+} (foam-like macrophage) or \code{fMAC-} (non-foam-like macrophage) and returns, for every cell, the \emph{FoamDEX}: the predicted probability of being a foam-like macrophage.
}

\usage{
FoamSpotteR(mac_obj, max.missing = 0.1)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{max.missing}{Maximum fraction of model genes that may be missing from the data. Missing genes are set to 0 with a warning; if more genes are missing the function stops. Default: \code{0.1}. Use \code{0} to require all the model genes.}
}

\details{
The organism-specific model (\code{\link{mouse_LAMprey}}, 20 genes, or \code{\link{human_LAMprey}}, 30 genes) is chosen from the \code{org} argument used in \code{\link{CreateMacObj}}. The model genes should be present in the expression matrix: avoid restricting the Seurat object to variable features before creating the MACanalyzeR object. Missing genes are treated as not expressed, which can lower the FoamDEX of foam-like cells.

Running \code{FoamSpotteR()} again overwrites the previous results.

The module is meant for myeloid cells, since the foam phenotype is typical of macrophages under specific conditions (e.g. lipid overload).
}

\value{
The input \code{MACanalyzeR} object with three new columns in \code{mac_obj@MetaData}:
  \item{fMAC+}{FoamDEX, the probability (between 0 and 1) that the cell is a foam-like macrophage.}
  \item{fMAC-}{Probability that the cell is not a foam-like macrophage (\code{1 - FoamDEX}).}
  \item{Foam}{Predicted class, \code{fMAC+} or \code{fMAC-}.}
}

\references{
Ninni A. \emph{et al.} (2025). MACanalyzeR scRNAseq analysis tool reveals PPARgamma(HIGH)/GDF15(HIGH) lipid-associated macrophages facilitate thermogenic expansion in BAT. \emph{Nature Communications} 16, 5063. \doi{10.1038/s41467-025-60295-2}
}

\seealso{\code{\link{FoamPlot}}, \code{\link{FoamLine}}, \code{\link{MacPlot}}}

\examples{
\dontrun{
mac_obj <- FoamSpotteR(mac_obj)
head(mac_obj@MetaData[, c("Foam", "fMAC+")])
table(mac_obj@MetaData$Foam)
}
}
