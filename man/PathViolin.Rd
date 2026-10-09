\name{PathViolin}
\alias{PathViolin}

\title{Violin plots of single-cell pathway scores}

\description{
Draws, for each selected pathway, a violin plot (with an overlaid boxplot) of the single-cell scores computed by \code{\link{PathAnalyzeR}} in each group, optionally with the p-values of all pairwise comparisons between groups.
}

\usage{
PathViolin(mac_obj, pathway, plot.by = mac_obj@ident, txt.size = 15, ncol = 1,
           col = NULL, intercept = FALSE, stat = FALSE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{PathAnalyzeR}} has been run with \code{meta = plot.by}.}
  \item{pathway}{Indices (as listed by \code{\link{PathDisplay}}) or names of one or more pathways. One violin plot is drawn for each pathway.}
  \item{plot.by}{Name of the grouping to compare. It must match the \code{meta} argument used in \code{\link{PathAnalyzeR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns used to arrange the plots when more than one pathway is given. Default: \code{1}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{intercept}{Logical. If \code{TRUE}, a horizontal reference line is drawn at y = 1 (average pathway activity). Default: \code{FALSE}.}
  \item{stat}{Logical. If \code{TRUE}, the p-values of the Wilcoxon rank-sum tests between all pairs of groups are added to the plots. Default: \code{FALSE}.}
}

\value{The plots are arranged in a grid and drawn with \code{\link[gridExtra]{grid.arrange}}, which also returns them invisibly as a \code{gtable} object.}

\seealso{\code{\link{PathAnalyzeR}}, \code{\link{PathDisplay}}}

\examples{
\dontrun{
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathDisplay(mac_obj, meta = "Sample")

PathViolin(mac_obj, pathway = c(1, 2), plot.by = "Sample", ncol = 2)
PathViolin(mac_obj, pathway = 1, plot.by = "Sample", intercept = TRUE, stat = TRUE)
}
}
