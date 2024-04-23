\name{unmasked-datasets}
\docType{data}
\alias{unmasked-datasets}
\alias{hg18Unmasked}
\alias{hg19Unmasked}
\alias{hg38Unmasked}
\alias{b36Unmasked}
\alias{b37Unmasked}
\title{Unmasked Regions of Human Genomes}
\description{Pre-built \code{\linkS4class{GRangesList}} objects with
  unmasked regions of different human genome builds}
\usage{hg18Unmasked
hg19Unmasked
hg38Unmasked
b36Unmasked
b37Unmasked}
\format{Each of these is a \code{\linkS4class{GRangesList}} object with
  unmasked regions of different human genome builds, as provided by the packages
  \pkg{BSgenome.Hsapiens.UCSC.hg18.masked},
  \pkg{BSgenome.Hsapiens.UCSC.hg19.masked}, and
  \pkg{BSgenome.Hsapiens.UCSC.hg38.masked}.
  The two latter,
  \code{b36Unmasked} and \code{b37Unmasked}, are variants using
  chromosome names as the genomes b36 and b37 that are frequently used
  by the Genome Analysis Toolkit (GATK).

  All four data sets comprise all 22 autosomal chromosomes, the two sex
  chromosomes, mitochondrial DNA, and the six pseudoautosomal regions as
  defined in the data frames
  \code{\link[GWASTools]{pseudoautosomal.hg18}} (for
  hg18), \code{\link[GWASTools]{pseudoautosomal.hg19}}
  (for hg19), and \code{\link[GWASTools]{pseudoautosomal.hg38}}
  (for hg38) as provided by the \pkg{GWASTools} package. If this is
  undesired, the user can re-unite the pseudoautosomal regions with
  their chromosomes as shown in the example section below or run
  \code{\link{unmaskedRegions}} him- or herself to extract unmasked regions.}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\link{unmaskedRegions}}, \code{\linkS4class{GRangesList}},
  \code{\linkS4class{GRanges}},
  \code{\link[GWASTools]{pseudoautosomal}}}
\examples{
## load data sets
data(hg19Unmasked)
data(b37Unmasked)

## show chromosome names
seqlevels(hg19Unmasked)
seqlevels(b37Unmasked)

## show names of list components
names(hg19Unmasked)
names(b37Unmasked)

## determine numbers of regions for each component
sapply(hg19Unmasked, length)
sapply(b37Unmasked, length)

## re-unite pseudoautosomal regions with their chromosomes
## (that is the result of what unmaskedRegions() would have given
## without specifying pseudoautosomal regions)
hg19basic <- hg19Unmasked[paste0("chr", 1:22)]
hg19basic$chrX <- reduce(unlist(hg19Unmasked[c("chrX", "X.PAR1",
                                               "X.PAR2", "X.XTR")]))
hg19basic$chrY <- reduce(unlist(hg19Unmasked[c("chrY", "Y.PAR1",
                                               "Y.PAR2", "Y.XTR")]))

## show some information about the newly created object
names(hg19basic)
sapply(hg19basic, length)
}
\keyword{datasets}

