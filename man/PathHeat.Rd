\name{PathHeat}
\alias{PathHeat}

\title{Heatmap of group-level pathway scores}

\description{
Draws a heatmap of the group-level pathway scores computed by \code{\link{PathAnalyzeR}}, giving an overview of which pathways are enriched in each group. Only pathways below the chosen p-value threshold are shown; rows are ordered by the group in which each pathway reaches its highest score.
}

\usage{
PathHeat(mac_obj, plot.by = mac_obj@ident, col = c("blue", "white", "red"),
         pval = 0.05, txt.size = 10, adjust = FALSE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{PathAnalyzeR}} has been run with \code{meta = plot.by}.}
  \item{plot.by}{Name of the grouping to plot. It must match the \code{meta} argument used in \code{\link{PathAnalyzeR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{col}{Vector of colours used to build the colour gradient, from low to high scores. Default: \code{c("blue", "white", "red")}.}
  \item{pval}{P-value threshold: only pathways with a p-value below it are plotted. Default: \code{0.05}.}
  \item{txt.size}{Font size. Default: \code{10}.}
  \item{adjust}{Logical. If \code{TRUE}, the threshold is applied to the Benjamini-Hochberg adjusted p-values (\code{padj}) instead of the raw ones. Default: \code{FALSE}.}
}

\details{
The colour scale is centred on 1 (average pathway activity) and is symmetric around it. Missing scores are shown as 1. If no pathway passes the threshold, a message is shown and nothing is plotted. With a single group (no p-values) all scored pathways are shown.
}

\value{The heatmap is drawn on the current graphics device and the \code{pheatmap} object is returned invisibly (see \code{\link[pheatmap]{pheatmap}}).}

\seealso{\code{\link{PathAnalyzeR}}, \code{\link{HeatMito}}}

\examples{
\dontrun{
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathHeat(mac_obj, plot.by = "Sample")

mac_obj <- PathAnalyzeR(mac_obj, meta = "Cluster")
PathHeat(mac_obj, plot.by = "Cluster", pval = 0.01)
}
}
