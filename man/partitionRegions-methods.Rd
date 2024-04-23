\name{partitionRegions-methods}
\docType{methods}
\alias{partitionRegions}
\alias{partitionRegions-methods}
\alias{partitionRegions,GRanges-method}
\alias{partitionRegions,GRangesList-method}
\alias{partitionRegions,MaskedBSgenome-method}
\title{Partition Genomic Regions}
\description{
  Partitions genomic regions into windows of fixed size
}
\usage{
\S4method{partitionRegions}{GRanges}(x, chrs=character(), width=5000, overlap=0.5)
\S4method{partitionRegions}{GRangesList}(x, chrs=character(), width=5000, overlap=0.5)
\S4method{partitionRegions}{MaskedBSgenome}(x, chrs=character(), width=5000, overlap=0.5, ...)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{GRanges}},
    \code{\linkS4class{GRangesList}}, or
    \code{\link[BSgenome:BSgenome]{MaskedBSgenome}}}
  \item{chrs}{a character vector (possibly empty) with names of
    chromosomes to limit to}
  \item{width}{window size}
  \item{overlap}{amount of overlap; a zero value corresponds to
    non-overlapping windows and the default 0.5 corresponds to 50\% overlap.
    The largest possible value is 0.8 which corresponds to an overlap of
    80\%.}
  \item{...}{further arguments are passed on to
    \code{\link{unmaskedRegions}}.}
}
\details{
  For a \code{\linkS4class{GRanges}} object \code{x}, this method
  partitions each genomic region into possibly overlapping, equally large
  windows of size \code{width}. The amount of overlap is controlled with
  the \code{overlap} parameter. The windows are placed such that
  possible overhangs are balanced at the beginning and end of the
  region. As an example, suppose we have a region from bases 1 to 14,000 and
  that we want to cover it with windows of 10,000bp length and
  50\% overlap. The straightforward approach would be to have two
  windows 1-10,000 and 5,001-15,000, and to crop the latter to
  5,001-14,000. As said, the \code{partitionRegions} balances the
  overhangs, so it will return two windows 1-9,500 and 4,501-14,000
  instead.
  
  If \code{chrs} is not empty, \code{partitionRegions} will only
  consider regions from those chromosomes (i.e. regions in the
  \code{\linkS4class{GRanges}} object whose \code{seqnames} occur
  in \code{chrs}). 

  If called for a \code{\linkS4class{GRangesList}} object, all
  componentes of the \code{\linkS4class{GRangesList}} object are
  partitioned separately as described above.

  For convenience, this function can also be called
  for a \code{\link[BSgenome:BSgenome]{MaskedBSgenome}} object. In this case,
  \code{\link{unmaskedRegions}} is called before partitioning.
}
\value{
  If \code{x} is a \code{\linkS4class{GRanges}} object, the function
  also returns a \code{\linkS4class{GRanges}} object. In the other two
  cases, a \code{\linkS4class{GRangesList}} object is returned.
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\link{assocTest}}, \code{\link{unmaskedRegions}},
  \code{\link{unmasked-datasets}}, \code{\linkS4class{GRangesList}},
  \code{\linkS4class{GRanges}}}
\examples{
## create a toy example
gr <- GRanges(seqnames="chr1", ranges=IRanges(start=1, end=14000))
partitionRegions(gr, width=10000, overlap=0.5)
partitionRegions(gr, width=10000, overlap=0.8)
partitionRegions(gr, width=10000, overlap=0)

## a toy example of a 'GRangesList'
grL <- GRangesList(gr, GRanges(seqnames="chr2",
                               ranges=IRanges(start=1, end=22000)))
partitionRegions(grL, width=10000, overlap=0.5)
partitionRegions(grL, width=10000, overlap=0.8)

## real-world example
data(hg38Unmasked)
partitionRegions(hg38Unmasked, chrs="chr22", width=20000)
}
\keyword{methods}
