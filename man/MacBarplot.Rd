\name{MacBarplot}
\alias{MacBarplot}

\title{Barplot of macrophage polarization classes}

\description{
Draws a stacked barplot with the percentage of \emph{Inflammatory}, \emph{Transitional} and \emph{Healing} cells (as assigned by \code{\link{MacPolarizeR}}) in each group.
}

\usage{
MacBarplot(mac_obj, plot.by = mac_obj@ident, txt.size = 14,
           col = c("#C5283D", "#E9724C", "#FFC857"))
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{MacPolarizeR}} has been run.}
  \item{plot.by}{Name of the metadata column defining the groups (one bar per group). Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{col}{Vector of three colours for the Inflammatory, Transitional and Healing classes, in this order. If a different number of colours is given, the default palette is used.}
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{MacPolarizeR}}, \code{\link{MacRadar}}}

\examples{
\dontrun{
mac_obj <- MacPolarizeR(mac_obj)

MacBarplot(mac_obj, plot.by = "Sample")
}
}
