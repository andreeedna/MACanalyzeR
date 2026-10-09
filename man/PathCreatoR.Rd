\name{PathCreatoR}
\alias{PathCreatoR}

\title{Retrieve MSigDB gene set collections}

\description{
Returns one of the MSigDB gene set collections bundled with the package, already formatted as a named list of gene vectors ready to be passed to the \code{pathway} argument of \code{\link{PathAnalyzeR}}. Gene sets can be filtered by name and size.
}

\usage{
PathCreatoR(organism, subcollection = NULL, pattern = NULL, min.size = 1, max.size = Inf)
}

\arguments{
  \item{organism}{Organism of the gene sets: \code{"mm"} (\emph{Mus musculus}) or \code{"hs"} (\emph{Homo sapiens}).}
  \item{subcollection}{Name of the MSigDB (sub)collection (see Details). Matching is case-insensitive, the last part of the name is enough when it is unique (e.g. \code{"GO:BP"} or \code{"KEGG_LEGACY"}), and human codes are accepted for mouse and vice versa (e.g. \code{"C5:GO:BP"} for \code{"M5:GO:BP"}). If \code{NULL}, the available collections are listed.}
  \item{pattern}{Optional regular expression: only gene sets whose name matches it (case-insensitive) are returned.}
  \item{min.size, max.size}{Minimum and maximum number of genes of the returned gene sets. Defaults: \code{1} and \code{Inf}.}
}

\details{
The gene sets are stored in the package (\code{\link{gsea_mouse}} and \code{\link{gsea_human}}), so no internet connection is needed. Available collections:

\itemize{
  \item Mouse (\code{"mm"}): \code{MH} (hallmark), \code{M1}, \code{M2:CGP}, \code{M2:CP}, \code{M2:CP:BIOCARTA}, \code{M2:CP:REACTOME}, \code{M2:CP:WIKIPATHWAYS}, \code{M3:GTRD}, \code{M3:MIRDB}, \code{M5:GO}, \code{M5:GO:BP}, \code{M5:GO:CC}, \code{M5:GO:MF}, \code{M5:MPT}, \code{M8}.
  \item Human (\code{"hs"}): \code{H} (hallmark), \code{C1}, \code{C2:CGP}, \code{C2:CP}, \code{C2:CP:BIOCARTA}, \code{C2:CP:REACTOME}, \code{C2:CP:WIKIPATHWAYS}, \code{C2:CP:PID}, \code{C2:CP:KEGG_LEGACY}, \code{C2:CP:KEGG_MEDICUS}, \code{C3:MIR}, \code{C3:MIR:MIRDB}, \code{C3:MIR:MIR_LEGACY}, \code{C3:TFT}, \code{C3:TFT:GTRD}, \code{C3:TFT:TFT_LEGACY}, \code{C4:3CA}, \code{C4:CGN}, \code{C4:CM}, \code{C5:GO}, \code{C5:GO:BP}, \code{C5:GO:CC}, \code{C5:GO:MF}, \code{C5:HPO}, \code{C6}, \code{C7:IMMUNESIGDB}, \code{C7:VAX}, \code{C8}.
}

Large collections (e.g. GO Biological Process) contain thousands of gene sets: use \code{pattern}, \code{min.size} and \code{max.size} to keep the analysis focused and fast.
}

\value{A named list of character vectors (one per gene set). When \code{subcollection} is \code{NULL}, the names of the available collections are returned invisibly.}

\references{
Liberzon A. \emph{et al.} (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. \emph{Cell Systems} 1(6), 417-425.
}

\seealso{\code{\link{PathAnalyzeR}}}

\examples{
\dontrun{
PathCreatoR("mm")                       # list the collections

hallmark <- PathCreatoR("mm", "MH")
macro_go <- PathCreatoR("mm", "GO:BP", pattern = "macrophage|foam", min.size = 10, max.size = 500)

mac_obj <- PathAnalyzeR(mac_obj, pathway = hallmark, meta = "Sample")
}
}
