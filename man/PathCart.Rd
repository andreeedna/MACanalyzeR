\name{PathCart}
\alias{PathCart}

\title{Cartplot of single-cell pathway scores}

\description{
The Cartplot is a visualization introduced by MACanalyzeR to look at the joint distribution of two pathways at the single-cell level. The single-cell scores of the two pathways computed by \code{\link{PathAnalyzeR}} are placed on the x and y axes and drawn as 2D density contours, one per group. Reference lines at average activity (0 on the default log2 scale, 1 on the original scale) split the plane into four quadrants, highlighting cells in which both pathways are co-upregulated. When a single pathway is given, a 1D density plot of its score is drawn instead.
}

\usage{
PathCart(mac_obj, pathway, plot.by = mac_obj@ident, split.by = NULL,
         fill = TRUE, a = 0.1, col = NULL, txt.size = 15, ncol = 3, scale = TRUE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{PathAnalyzeR}} has been run with \code{meta = plot.by}.}
  \item{pathway}{Indices (as listed by \code{\link{PathDisplay}}) or names of one or two pathways.}
  \item{plot.by}{Name of the grouping used to colour the densities. It must match the \code{meta} argument used in \code{\link{PathAnalyzeR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{fill}{Logical. If \code{TRUE} the density contours are filled, otherwise only the contour lines are drawn. Default: \code{TRUE}.}
  \item{a}{Transparency (alpha) of the filled densities, between 0 and 1. Default: \code{0.1}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{scale}{Logical. If \code{TRUE}, scores are shown on a log2 scale and the axis labels become \code{log2(...)}: the reference line moves from 1 to 0, and values of +1 and -1 mean twice and half the average activity. Cells with a score of 0 or below cannot be log-transformed and are removed (their number is reported). If \code{FALSE}, the original scores are shown. Default: \code{TRUE}.}
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{PathAnalyzeR}}, \code{\link{PathDisplay}}, \code{\link{PolCart}}}

\examples{
\dontrun{
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathDisplay(mac_obj, meta = "Sample")

## two pathways: 2D Cartplot
PathCart(mac_obj, pathway = c(1, 2), plot.by = "Sample")

## one pathway: 1D density
PathCart(mac_obj, pathway = 1, plot.by = "Sample", fill = FALSE)

## original (linear) scores
PathCart(mac_obj, pathway = c(1, 2), plot.by = "Sample", scale = FALSE)
}
}
