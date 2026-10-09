\name{macSpectrum}
\alias{macSpectrum}
\alias{MacSpectrumIndex}

\title{MacSpectrum and AtheroSpectrum plots}

\description{
Computes, for every cell, macrophage activation indices based on the correlation of its transcriptome with reference macrophage profiles, and plots the cells on the resulting two-dimensional spectrum. Two modes are available:
\itemize{
  \item \code{"macspectrum"}: macrophage polarization index (MPI, x axis) versus activation-induced macrophage differentiation index (AMDI, y axis), following the MacSpectrum approach;
  \item \code{"atherospectrum"}: MPI (x axis) versus a macrophage foaming differentiation index (MDFI, y axis), which measures how closely each cell resembles foamy rather than control macrophages.
}
}

\usage{
macSpectrum(mac_obj, mode = "macspectrum", plot.by = mac_obj@ident, split.by = NULL,
            style = "density", txt.size = 17, col = NULL, ncol = 3)

MacSpectrumIndex(mac_obj)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{mode}{Type of spectrum to plot: \code{"macspectrum"} or \code{"atherospectrum"}. Default: \code{"macspectrum"}.}
  \item{plot.by}{Name of the metadata column used to colour the cells. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{style}{\code{"density"} to draw filled 2D density contours or \code{"point"} to draw one point per cell. Default: \code{"density"}.}
  \item{txt.size}{Base text size. Default: \code{17}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. If its length does not match the number of groups a warning is raised and the default palette is used.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
}

\details{
\itemize{
  \item \strong{MPI} places each cell between the M1 (positive values) and M2 (negative values) reference profiles, on a scale from -50 to 50.
  \item \strong{AMDI} is derived from the (negative) correlation with the reference M0 profile: higher values indicate cells that are more differentiated/activated.
  \item \strong{MDFI} places each cell between control (negative values) and foamy (positive values) reference macrophages, on a scale from -50 to 50.
}
The reference profiles (\code{\link{M0_mean_new}}, \code{\link{M1_mean_new}}, \code{\link{M2_mean_new}}, \code{\link{con_mean_new}}, \code{\link{foam_mean_new}}) are defined on mouse genes; for human data genes are converted with \code{\link{human_mouse_dict}}.

\code{macSpectrum()} computes the indices on the fly, unless they are already stored in the object. \code{MacSpectrumIndex()} computes them and stores them in \code{mac_obj@MetaData} (columns \code{MPI}, \code{AMDI} and \code{MDFI}), so they can be reused, exported or plotted with other functions.
}

\value{\code{macSpectrum()}: a \code{ggplot} object. \code{MacSpectrumIndex()}: the input object with the indices added to \code{mac_obj@MetaData}.}

\references{
Li C. \emph{et al.} (2019). Single-cell transcriptomics-based MacSpectrum reveals novel macrophage activation signatures in diseases. \emph{JCI Insight} 4(10), e126453.
}

\seealso{\code{\link{MacPolarizeR}}, \code{\link{FoamSpotteR}}}

\examples{
\dontrun{
macSpectrum(mac_obj, plot.by = "Sample")
macSpectrum(mac_obj, mode = "atherospectrum", plot.by = "Cluster", style = "point")

mac_obj <- MacSpectrumIndex(mac_obj)
summary(mac_obj@MetaData[, c("MPI", "AMDI", "MDFI")])
}
}
