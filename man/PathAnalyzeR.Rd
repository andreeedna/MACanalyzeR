\name{PathAnalyzeR}
\alias{PathAnalyzeR}

\title{
Pathway Enrichment Analysis
}
\description{
The PathAnalyzeR module enables the enrichment analysis of cells based on gene expression. The function PathAnalyzeR() performs an enrichment analysis for each cell individually and computes a score for each Cluster or Sample. The technique for calculating the enrichment score for each gene set is elaborated in the MACanalyzeR paper.
}
\usage{
mac_obj <- PathAnalyzeR(mac_obj,
                        pathway = NULL,
                        meta = "Sample"
                       )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{pathway}{
a list of pathways of interest. Each element must be a vector containing the genes involved in a pathway. These vectors must be named with the corresponding pathway names. By default, it uses a subset of KEGG pathways of interest in biological processes or cellular components of macrophages.
}
  \item{meta}{
metadata selected to perform PathAnalyzeR
}
}
\details{
A collection of KEGG, GO and Reactome for mouse (*Mus musculus*) is already incorporated into the package. However, users have the flexibility to utilize any category extracted from other databases by providing a list of genes with a named pathway associated.
}

\note{
PathAnalyzeR allows you to define custom pathways for analysis. These pathways can be:

- Subsets of existing databases like KEGG, Gene Ontology (GO), and Reactome pathways (which are already implemented as lists in MACanalyzeR);

- User-created pathways: you can build your own pathway by specifying a list of genes;

}

\examples{
mac_obj <- PathAnalyzeR(mac_obj)

# subset of KEGG, GO or REACTOME
ptw <- kegg[c(1:10)]
ptw <- go[c(1:10)]
ptw <- reactome[c(1:10)]

# user-created pathways
ptw <- list("path_name1" = c(...),
			"path_name2" = c(...),
			...
			)

# PathAnalyzeR performed on Cluster
mac_obj <- PathAnalyzeR(mac_obj, pathway=ptw, meta="Cluster")
}
