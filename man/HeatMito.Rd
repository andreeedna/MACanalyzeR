\name{HeatMito}
\alias{HeatMito}

\title{Heatmaps of mitochondrial pathway scores}

\description{
Draws the results of \code{\link{MitoScanneR}} as a series of heatmaps, one for each MitoCarta category, showing the score of every mitochondrial pathway in each group. Rows are ordered by the group in which each pathway reaches its highest score.
}

\usage{
HeatMito(mac_obj, plot.by = mac_obj@ident, col = c("blue", "white", "red"),
         txt.size = 10)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{MitoScanneR}} has been run with \code{meta = plot.by}.}
  \item{plot.by}{Name of the grouping to plot. It must match the \code{meta} argument used in \code{\link{MitoScanneR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{col}{Vector of colours used to build the colour gradient, from low to high scores. Default: \code{c("blue", "white", "red")}.}
  \item{txt.size}{Font size. Default: \code{10}.}
}

\details{
The colour scale of each heatmap is centred on 1 (average activity) and is symmetric around it. Pathways that could not be scored (\code{NA}) are shown as 1.
}

\value{The heatmaps are drawn on the current graphics device, one after the other (see \code{\link[pheatmap]{pheatmap}}).}

\seealso{\code{\link{MitoScanneR}}, \code{\link{PathHeat}}}

\examples{
\dontrun{
mac_obj <- MitoScanneR(mac_obj, meta = "Sample")
HeatMito(mac_obj, plot.by = "Sample")
}
}
