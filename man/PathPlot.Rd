\name{PathPlot}
\alias{PathPlot}

\title{Plot a single-cell pathway score on the dimensionality reduction}

\description{
Draws a 2D scatter plot of the cells on the embedding, coloured by their single-cell score for one or more pathways, as computed by \code{\link{PathAnalyzeR}}. This shows where a pathway is active across the cell landscape. Cells with the highest scores are drawn on top.
}

\usage{
PathPlot(mac_obj, pathway, txt.size = 14, shade = "turbo", split.by = NULL,
         pt.size = 2, ncol = 2, max.cutoff = NA, meta = mac_obj@ident)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{PathAnalyzeR}} has been run with \code{meta = meta}.}
  \item{pathway}{Indices (as listed by \code{\link{PathDisplay}}) or names of one or more pathways. With more than one pathway, each is drawn in its own panel.}
  \item{txt.size}{Base text size. Default: \code{14}.}
  \item{shade}{Viridis colour scale used for the score. One of \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"viridis"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"} or \code{"turbo"}. Default: \code{"turbo"}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{pt.size}{Point size. Default: \code{2}.}
  \item{ncol}{Number of columns of the panel grid. Default: \code{2}.}
  \item{max.cutoff}{Optional upper limit of the score: higher values are set to this value, so that a few outlier cells do not compress the colour scale. Values below 2 are not recommended. Default: \code{NA} (no cutoff).}
  \item{meta}{Name of the grouping whose PathAnalyzeR results are used. It must match the \code{meta} argument used in \code{\link{PathAnalyzeR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
}

\details{
Single-cell scores depend on the grouping used by \code{\link{PathAnalyzeR}}, because each gene is compared with its average across groups. With several pathways and \code{split.by}, pathways are shown in rows and the groups of \code{split.by} in columns.
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{PathAnalyzeR}}, \code{\link{PathDisplay}}, \code{\link{MacPlot}}}

\examples{
\dontrun{
mac_obj <- PathAnalyzeR(mac_obj)
PathDisplay(mac_obj)

PathPlot(mac_obj, pathway = 1)
PathPlot(mac_obj, pathway = 1, split.by = "Sample", max.cutoff = 3)
PathPlot(mac_obj, pathway = c(1, 2))

mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathPlot(mac_obj, pathway = "OXIDATIVE_PHOSPHORYLATION", meta = "Sample")
}
}
