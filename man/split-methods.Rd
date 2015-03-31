\name{split-methods}
\docType{methods}
\alias{split}
\alias{split-methods}
\alias{split,GRanges,GRangesList-method}
\title{Split \code{GRanges} Object}
\description{
  Splits a \code{\linkS4class{GRanges}} object into a
  \code{\linkS4class{GRangesList}}
}
\usage{
\S4method{split}{GRanges,GRangesList}(x, f)
}
\arguments{
  \item{x}{object of class \code{\linkS4class{GRanges}}}
  \item{f}{object of class \code{\linkS4class{GRangesList}}}
}
\details{This function splits a \code{\linkS4class{GRanges}} object
  \code{x} along a \code{\linkS4class{GRangesList}} object \code{f}.
  More specifically, each region in \code{x} is checked for overlaps
  with every list component of \code{f}. The function returns a
  \code{\linkS4class{GRangesList}} object each component of which
  contains all overlaps of \code{x} with one of the components of
  \code{f}. If the overlap is empty, this component is discarded.

  This function is mainly made for splitting regions of interests
  (transcripts, exons, regions targeted by exome capturing)
  along chromosomes (and pseudoautosomal regions).

  The returned object inherits sequence infos (chromosome names,
  chromosome lengths, genome, etc.) from the
  \code{\linkS4class{GRangesList}} object \code{f}.

  For greater universality, the function takes strand information
  into account. If overlaps should not be determined in a
  strand-specific manner, all strand information must be discarded from
  \code{x} and \code{f} before calling \code{split}.
}
\value{
  a \code{\linkS4class{GRangesList}} object (see details above)
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\linkS4class{GRanges}}, \code{\linkS4class{GRangesList}}}
\examples{
## set up toy example
chr1 <- GRanges(seqnames="chr1", ranges=IRanges(start=1, end=200000))
chr2 <- GRanges(seqnames="chr2", ranges=IRanges(start=1, end=180000))
grL <- GRangesList(list(chr1=chr1, chr2=chr2))
seqlevels(grL) <- c("chr1", "chr2")
seqlengths(grL) <- c(chr1=200000, chr2=180000)
grL

## split set of regions given as 'GRanges' object
gr <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr2", "chr2"),
              ranges=IRanges(start=c(1, 30000, 10000, 51000, 110000),
                             end=c(340, 37000, 10100, 61000, 176000)))
gr
split(gr, grL)

## consider transcripts on the X chromosome, but with pseudoautosomal
## regions treated separately
if (require(TxDb.Hsapiens.UCSC.hg38.knownGene))
{
    data(hg38Unmasked)

    hg38tr <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    strand(hg38tr) <- "*"

    split(hg38tr, hg38Unmasked[c("chrX", "X.PAR1", "X.PAR2", "X.XTR")])
}}
\keyword{methods}
