\name{PathDisplay}
\alias{PathDisplay}

\title{List the pathways analysed by PathAnalyzeR}

\description{
Prints the pathways analysed by \code{\link{PathAnalyzeR}} for a given grouping, each preceded by its index number. The single-cell plotting functions (\code{\link{PathCart}}, \code{\link{PathPlot}} and \code{\link{PathViolin}}) select pathways by these numbers.
}

\usage{
PathDisplay(mac_obj, meta = mac_obj@ident)
}

\arguments{
  \item{mac_obj}{A \code{MACanalyzeR} object on which \code{\link{PathAnalyzeR}} has been run.}
  \item{meta}{Name of the grouping whose results should be listed. It must match the \code{meta} argument used in \code{\link{PathAnalyzeR}}. Default: the grouping set with the \code{id} argument of \code{\link{CreateMacObj}}.}
}

\value{The numbered list of pathways is printed as messages. A data frame with the columns \code{index} and \code{pathway} is returned invisibly.}

\seealso{\code{\link{PathAnalyzeR}}, \code{\link{PathCart}}, \code{\link{PathPlot}}, \code{\link{PathViolin}}}

\examples{
\dontrun{
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathDisplay(mac_obj, meta = "Sample")
}
}
