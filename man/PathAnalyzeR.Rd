\name{PathAnalyzeR}
\alias{PathAnalyzeR}

\title{Pathway enrichment scores for groups and single cells}

\description{
The PathAnalyzeR module scores the activity of a set of pathways (gene sets) both at the group level (e.g. per sample or cluster) and in each single cell, allowing pathway activity to be compared across conditions and visualized on the embedding.
}

\usage{
PathAnalyzeR(mac_obj, pathway = NULL, meta = NULL, verbose = TRUE)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object created with \code{\link{CreateMacObj}}.}
  \item{pathway}{A named list of gene sets: each element is a character vector of gene symbols and its name is the pathway name. A single character vector is also accepted and is analysed as one gene set called \code{"Custom gene set"}. Gene symbols must match the organism of the dataset. Default: \code{NULL}, which uses a curated selection of KEGG pathways relevant to macrophage biology (\code{\link{mouse_selected}} or \code{\link{human_selected}}).}
  \item{meta}{Name of the metadata column defining the groups to compare. Default: \code{NULL}, which uses the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
  \item{verbose}{Logical. If \code{TRUE}, a progress bar and messages are shown. Default: \code{TRUE}.}
}

\details{
For each pathway, the expression of every gene is divided by its average expression across the groups of \code{meta}; the resulting ratios are then averaged, weighting each gene by its expression variance. A score of 1 means average activity, values above 1 indicate higher and values below 1 lower activity than the dataset average. The full method is described in the MACanalyzeR paper.

Pathways are skipped (scores set to \code{NA}) when fewer than 5 of their genes are detected, or when fewer than 3 genes have a mean expression of at least 0.001 in every group.

Ready-to-use gene set collections are provided with the package: \code{\link{mouse_kegg}}, \code{\link{mouse_go}}, \code{\link{mouse_reactome}}, and the MSigDB collections for mouse and human available through \code{\link{PathCreatoR}}. Any custom gene set can be used as well.

The p-value tests whether the genes of a pathway are consistently shifted between groups. With two groups, it is a one-sample Wilcoxon signed-rank test of the gene ratios of the first group against 1 (equivalent to comparing the two groups). With three or more groups, it is a Friedman test with genes as blocks. Before version 1.1.0 the two-group test was also applied, to the first group only, when more groups were present; p-values for three or more groups therefore differ from earlier versions. With a single group no p-value is computed.

Results are stored under the name of \code{meta}, so the function can be run several times with different groupings and all results are kept. Plotting functions look up the results by the same name (e.g. \code{PathHeat(mac_obj, plot.by = "Sample")} requires \code{PathAnalyzeR(mac_obj, meta = "Sample")}).
}

\value{
The input \code{MACanalyzeR} object with a new element \code{mac_obj@PathAnalyzeR[[meta]]}, a list of two data frames:
  \item{Total}{Pathways x groups matrix of group-level scores, plus a \code{pval} column testing whether the pathway differs between groups (see Details) and a \code{padj} column with the Benjamini-Hochberg adjusted p-values.}
  \item{SingleCell}{Cells x pathways matrix of single-cell scores.}
}

\references{
Ninni A. \emph{et al.} (2025). MACanalyzeR scRNAseq analysis tool reveals PPARgamma(HIGH)/GDF15(HIGH) lipid-associated macrophages facilitate thermogenic expansion in BAT. \emph{Nature Communications} 16, 5063. \doi{10.1038/s41467-025-60295-2}
}

\seealso{
\code{\link{PathHeat}}, \code{\link{PathDisplay}}, \code{\link{PathCart}}, \code{\link{PathPlot}}, \code{\link{PathViolin}}, \code{\link{PathCreatoR}}
}

\examples{
\dontrun{
## default macrophage-relevant KEGG pathways, grouped by sample
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")

## a subset of a built-in collection
ptw <- mouse_kegg[1:10]

## an MSigDB collection
ptw <- PathCreatoR("mm", "MH")

## custom gene sets
ptw <- list("Lipid uptake"  = c("Cd36", "Lpl", "Msr1", "Scarb1", "Ldlr"),
            "LAM signature" = c("Trem2", "Gpnmb", "Lgals3", "Cd9", "Spp1"))

mac_obj <- PathAnalyzeR(mac_obj, pathway = ptw, meta = "Cluster")
head(mac_obj@PathAnalyzeR$Cluster$Total)
}
}
