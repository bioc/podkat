\name{hgA}
\docType{data}
\alias{hgA}
\title{Artificial Human Chromosome for Testing Purposes}
\description{A \code{\linkS4class{GRanges}} object defining a
  minimalistic artificial human chromosome with 200,000 bp length}
\usage{hgA}
\format{Real human genome-based examples would require the supply
  of massive data and would require lengthy computation
  times. Therefore, the examples supplied with this package are based on
  a small single-chromosome artificial genome. The
  \code{\linkS4class{GRanges}} object \code{hgA} provides a description
  of this artificial genome that can be used for further processing,
  e.g. by the \code{\link{partitionRegions}} function.}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\linkS4class{GRanges}}, \code{\link{partitionRegions}}}
\examples{
## load data set
data(hgA)

## display hgA
show(hgA)
genome(hgA)

## partition hgA into overlapping regions of length 10,000 bp
partitionRegions(hgA, width=10000)
}
\keyword{datasets}

