\name{MacPolarizeR}
\alias{MacPolarizeR}

\title{Classify macrophage polarization}

\description{
The MacPolarizeR module classifies monocytes and macrophages according to their polarization state. Each cell receives an M1 (pro-inflammatory) and an M2 (anti-inflammatory) score computed from curated marker genes; cells are then grouped into three states along the M1/M2 polarization balance: \emph{Inflammatory} (M1-like), \emph{Transitional} (intermediate) and \emph{Healing} (M2-like).
}

\usage{
MacPolarizeR(mac_obj, seed = 42, nstart = 25, n.boot = 20)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{seed}{Random seed used for the k-means initialization and for the bootstrap, so that results are reproducible. The random number generator state of the session is left unchanged. Use \code{NULL} to use the current state instead. Default: \code{42}.}
  \item{nstart}{Number of random starts of k-means; the best solution is kept. Default: \code{25}.}
  \item{n.boot}{Number of bootstrap replicates used to estimate the stability of each assignment (see Details). Use \code{0} to skip it. Default: \code{20}.}
}

\details{
The algorithm works in four steps:
\enumerate{
  \item \strong{Scores.} The M1 and M2 scores of each cell are computed from the marker genes (\code{\link{mouse_PolGenes}} or \code{\link{human_PolGenes}}, chosen from the \code{org} argument of \code{\link{CreateMacObj}}) with the MACanalyzeR scoring method (see \code{\link{PathAnalyzeR}}), relative to the average of the whole dataset: a score of 1 means average expression of the gene set. Scores do not depend on any clustering.
  \item \strong{Polarization balance.} Both scores are standardized and their difference, \code{z(M1) - z(M2)}, gives the polarization balance of each cell. Using the difference removes the component shared by the two scores (overall activation, transcriptional activity or cell quality), which would otherwise dominate the classification.
  \item \strong{Classification.} The balance is split into three groups by k-means with \code{nstart} random starts. Being one-dimensional, the solution is stable across runs. Groups are named by their centroid: highest balance \emph{Inflammatory}, lowest \emph{Healing}, intermediate \emph{Transitional}.
  \item \strong{Stability.} Marker genes are resampled with replacement \code{n.boot} times and the classification is repeated. The stability of a cell is the fraction of replicates assigning it to the same state; low values flag cells close to the boundary between two states.
}

The three states are defined relative to the dataset: they describe the polarization spectrum of the cells being analysed, and the same cell may be classified differently in a dataset with a different composition.

Version 1.1.0 replaced the previous implementation, which clustered cells by k-means directly on the marker gene expression with a single random start. Results from earlier versions are therefore not identical and were not reproducible across runs.
}

\value{
The input \code{MACanalyzeR} object with:
  \item{MetaData$Mac}{A factor with the polarization state of each cell (\code{Inflammatory}, \code{Transitional}, \code{Healing}).}
  \item{MetaData$Mac.stability}{Fraction of bootstrap replicates agreeing with the assignment (only when \code{n.boot > 0}).}
  \item{MacPolarizeR$SingleCell}{A matrix (cells x 2) with the M1 and M2 scores of each cell.}
  \item{MacPolarizeR$Total}{A matrix (2 x 3) with the mean M1 and M2 scores of each state.}
  \item{MacPolarizeR$Balance}{The polarization balance of each cell.}
  \item{MacPolarizeR$Genes}{The M1 and M2 genes used for the scores.}
  \item{MacPolarizeR$Params}{The parameters used, for reproducibility.}
}

\seealso{\code{\link{MacPlot}}, \code{\link{MacBarplot}}, \code{\link{MacRadar}}, \code{\link{PolCart}}}

\examples{
\dontrun{
mac_obj <- MacPolarizeR(mac_obj)

table(mac_obj@MetaData$Mac)
summary(mac_obj@MetaData$Mac.stability)
head(mac_obj@MacPolarizeR$SingleCell)

## keep only confidently assigned cells
confident <- mac_obj@MetaData$Mac.stability >= 0.8
}
}
