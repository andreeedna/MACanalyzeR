\name{GeneBalance}
\alias{GeneBalance}

\title{Balance between two gene sets in single cells}

\description{
Compares, in each cell, the activity of two gene sets and shows the distribution of their balance in each group as a violin plot. It is the general function behind \code{\link{MitoBalance}} and \code{\link{GliOxBalance}}, and can be used with any pair of gene sets (e.g. glycolysis vs fatty acid oxidation, lipid uptake vs efflux).
}

\usage{
GeneBalance(mac_obj, set1, set2, plot.by = mac_obj@ident, split.by = NULL, txt.size = 15,
            ncol = 3, col = NULL, intercept = FALSE, title = NULL,
            ylab = "log10(set1 / set2) score ratio")
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{set1, set2}{Character vectors with the genes of the two gene sets.}
  \item{plot.by}{Name of the metadata column defining the groups to compare. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{intercept}{Logical. If \code{TRUE}, a horizontal reference line is drawn at y = 0, where the two gene sets are balanced. Default: \code{FALSE}.}
  \item{title}{Optional plot title. Default: \code{NULL}.}
  \item{ylab}{Label of the y axis.}
}

\details{
The single-cell scores of the two gene sets are computed as in \code{\link{PathAnalyzeR}}, using the groups of \code{plot.by}. The value plotted for each cell is \code{log10(score set1 / score set2)}: positive values indicate a relative prevalence of \code{set1}, negative values of \code{set2}. Cells with non-finite values (a score equal to 0) are removed, and their number is reported.

Each gene set needs at least 5 genes detected in the data, of which at least 3 are expressed in every group; otherwise the function stops with an error reporting how many genes were found.
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{MitoBalance}}, \code{\link{GliOxBalance}}, \code{\link{PathAnalyzeR}}}

\examples{
\dontrun{
sel <- MACanalyzeR::mouse_selected
GeneBalance(mac_obj, sel$GLYCOLYSIS_GLUCONEOGENESIS, sel$FATTY_ACID_METABOLISM,
            plot.by = "Sample", intercept = TRUE, ylab = "log10(glycolysis / FAO)")
}
}
