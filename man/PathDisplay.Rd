\name{PathDisplay}
\alias{PathDisplay}

\title{
Pathway Display
}

\description{
In PathAnalyzeR enrichment single-cell visualization, all functions use pathwat identifiers in form of numbers. To estabilish this association between pathways and their corrisponding numbers, the PathDisplay() function is implemented.
}

\usage{
PathDisplay(mac_obj, meta = "Sample")
}

\arguments{
  \item{mac_obj}{
MACanalyzeR object
}
  \item{meta}{
metadata selected to perform PathDisplay
}
}

\examples{
PathDisplay(mac_obj)
}
