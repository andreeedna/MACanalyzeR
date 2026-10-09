\name{MitoBalance}
\alias{MitoBalance}

\title{Mito-nuclear balance of OXPHOS genes}

\description{
Compares, in each cell, the activity of the OXPHOS subunits encoded by the mitochondrial genome with that of the OXPHOS subunits encoded by the nuclear genome, and shows the distribution of this mito-nuclear balance in each group as a violin plot. Shifts in the balance can reveal mito-nuclear imbalance or altered mitochondrial biogenesis.
}

\usage{
MitoBalance(mac_obj, plot.by = mac_obj@ident, txt.size = 15, ncol = 3, col = NULL,
            intercept = FALSE, title = "Mito-nuclear balance", split.by = NULL)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{plot.by}{Name of the metadata column defining the groups to compare. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{txt.size}{Base text size. Default: \code{15}.}
  \item{ncol}{Number of columns of the panel grid when \code{split.by} is used. Default: \code{3}.}
  \item{col}{Optional vector of colours, one per group of \code{plot.by}. It must have exactly one colour per group, otherwise an error is raised.}
  \item{intercept}{Logical. If \code{TRUE}, a horizontal reference line is drawn at y = 0, where mitochondrial and nuclear subunits are balanced. Default: \code{FALSE}.}
  \item{title}{Plot title. Default: \code{"Mito-nuclear balance"}.}
  \item{split.by}{Optional name of a metadata column used to split the plot into panels. Default: \code{NULL} (no split).}
}

\details{
The MitoCarta 3.0 gene sets \emph{OXPHOS mitochondrial subunits} and \emph{OXPHOS nuclear subunits} are compared with \code{\link{GeneBalance}}: the value plotted for each cell is \code{log10(mitochondrial score / nuclear score)}. Positive values indicate a relative excess of mitochondrially encoded subunits, negative values of nuclear-encoded ones.

The 13 mitochondrially encoded genes (\code{mt-} genes in mouse, \code{MT-} in human) are often removed from single-cell data during quality control. They must be kept in the Seurat object to use this function, which otherwise stops with an informative error.

This function does not require \code{\link{MitoScanneR}} to be run first.
}

\value{A \code{ggplot} object.}

\seealso{\code{\link{GeneBalance}}, \code{\link{GliOxBalance}}, \code{\link{MitoScanneR}}}

\examples{
\dontrun{
MitoBalance(mac_obj, plot.by = "Sample", intercept = TRUE)
}
}
