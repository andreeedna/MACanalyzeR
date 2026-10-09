\name{MitoScanneR}
\alias{MitoScanneR}

\title{Mitochondrial and oxidative profiling}

\description{
The MitoScanneR module characterizes the mitochondrial and oxidative profile of the cells. It applies the PathAnalyzeR scoring approach to the pathways of MitoCarta 3.0, the inventory of mammalian mitochondrial proteins and pathways, and scores each mitochondrial pathway in each group, highlighting groups with higher or lower mitochondrial activity.
}

\usage{
MitoScanneR(mac_obj, meta = mac_obj@ident, verbose = TRUE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{meta}{Name of the metadata column defining the groups to compare. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{verbose}{Logical. If \code{TRUE}, progress bars and messages are shown. Default: \code{TRUE}.}
}

\details{
The organism-specific MitoCarta pathways (\code{\link{mouse_MitoCarta}} or \code{\link{human_MitoCarta}}) are chosen from the \code{org} argument used in \code{\link{CreateMacObj}}. Pathways are organized in the main MitoCarta categories (e.g. OXPHOS, Metabolism, Mitochondrial central dogma, Protein import, sorting and homeostasis), and each category is scored separately.

Scores are computed as in \code{\link{PathAnalyzeR}}: a value of 1 means average activity, values above 1 indicate higher activity. Pathways with fewer than 5 detected genes, or fewer than 3 genes expressed in every group, are set to \code{NA}.
}

\value{
The input \code{MACanalyzeR} object with a new element \code{mac_obj@MitoScanneR[[meta]]}: a list with one data frame (pathways x groups) per MitoCarta category.
}

\references{
Rath S. \emph{et al.} (2021). MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. \emph{Nucleic Acids Research} 49(D1), D1541-D1547.
}

\seealso{\code{\link{HeatMito}}, \code{\link{MitoBalance}}, \code{\link{GliOxBalance}}, \code{\link{PathAnalyzeR}}}

\examples{
\dontrun{
mac_obj <- MitoScanneR(mac_obj, meta = "Sample")
names(mac_obj@MitoScanneR$Sample)
HeatMito(mac_obj, plot.by = "Sample")
}
}
