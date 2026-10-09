\name{GliOxBalance}
\alias{GliOxBalance}

\title{Glycolysis/OXPHOS balance}

\description{
Compares, in each cell, the activity of glycolysis with that of oxidative phosphorylation (OXPHOS) and shows the distribution of this metabolic balance in each group as a violin plot. It helps to detect shifts towards a glycolytic (typically pro-inflammatory) or oxidative (typically pro-resolving) metabolism.
}

\usage{
GliOxBalance(mac_obj, plot.by = mac_obj@ident, txt.size = 15, ncol = 3, col = NULL,
             intercept = FALSE, title = "Glycolysis/OXPHOS balance", split.by = NULL)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{plot.by}{Name of the metadata column defining the groups to compare. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{intercept}{Logical. If \code{TRUE}, a horizontal reference line is drawn at y = 0, where glycolysis and OXPHOS are balanced. Default: \code{FALSE}.}
  \item{title}{Plot title. Default: \code{"Glycolysis/OXPHOS balance"}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
}

\details{
The KEGG pathways \emph{Glycolysis / Gluconeogenesis} and \emph{Oxidative phosphorylation} (from \code{\link{mouse_selected}} or \code{\link{human_selected}}) are compared with \code{\link{GeneBalance}}: the value plotted for each cell is \code{log10(glycolysis score / OXPHOS score)}. Positive values indicate a relatively glycolytic profile, negative values a relatively oxidative one.

This function does not require \code{\link{MitoScanneR}} to be run first.
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{GeneBalance}}, \code{\link{MitoBalance}}, \code{\link{MitoScanneR}}}

\examples{
\dontrun{
GliOxBalance(mac_obj, plot.by = "Sample", intercept = TRUE)
GliOxBalance(mac_obj, plot.by = "Mac", split.by = "Sample")
}
}
