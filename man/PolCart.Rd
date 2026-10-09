\name{PolCart}
\alias{PolCart}

\title{Cartplot of M1 and M2 polarization scores}

\description{
Shows the single-cell M1 and M2 polarization scores computed by \code{\link{MacPolarizeR}} on a Cartesian plane (M1 score on the x axis, M2 score on the y axis), either as 2D density contours or as points, coloured by group. Reference lines at average activity (0 on the default log2 scale, 1 on the original scale) split the plane into four quadrants (low/high M1 and low/high M2).
}

\usage{
PolCart(mac_obj, plot.by = mac_obj@ident, split.by = NULL, style = "density",
        fill = TRUE, a = 0.1, col = NULL, txt.size = 15, ncol = 3, scale = TRUE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{MacPolarizeR}} has been run.}
  \item{plot.by}{Name of the metadata column used to colour the cells. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{style}{\code{"density"} to draw filled 2D density contours or \code{"point"} to draw one point per cell. Default: \code{"density"}.}
  \item{fill}{Logical. With \code{style = "density"}, if \code{TRUE} the density contours are filled, otherwise only the contour lines are drawn. Default: \code{TRUE}.}
  \item{a}{Transparency (alpha) of the filled density contours, between 0 and 1. Default: \code{0.1}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{scale}{Logical. If \code{TRUE}, scores are shown on a log2 scale and the axis labels become \code{log2(...)}: the reference line moves from 1 to 0, and values of +1 and -1 mean twice and half the average activity. Cells with a score of 0 or below cannot be log-transformed and are removed (their number is reported). If \code{FALSE}, the original scores are shown. Default: \code{TRUE}.}
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{MacPolarizeR}}, \code{\link{PathCart}}}

\examples{
\dontrun{
mac_obj <- MacPolarizeR(mac_obj)

PolCart(mac_obj, plot.by = "Sample")
PolCart(mac_obj, plot.by = "Mac", style = "point")
PolCart(mac_obj, plot.by = "Sample", scale = FALSE)   # original scores
}
}
