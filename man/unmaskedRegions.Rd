\name{unmaskedRegions}
\alias{unmaskedRegions}
\title{Extract Unmasked Regions from \code{MaskedBSgenome} Object}
\description{
  Create a \code{\linkS4class{GRangesList}} of unmasked regions from
  a \code{\link[BSgenome:BSgenome]{MaskedBSgenome}} object
}
\usage{
unmaskedRegions(x, chrs=character(), pseudoautosomal=NULL,
                ignoreGaps=250, activeMasks=active(masks(x[[1]])))
}
\arguments{
  \item{x}{a \code{\link[BSgenome:BSgenome]{MaskedBSgenome}} object}
  \item{chrs}{a character vector of chromosome names to restrict to; if
     empty (default), all chromosomes in \code{x} are considered.}
  \item{pseudoautosomal}{if \code{NULL} (default), the chromosomes are
   considered as they are; \code{pseudoautosomal} must be a data
   frame complying with the format of the
   \code{\link[GWASTools]{pseudoautosomal.hg18}},
   \code{\link[GWASTools]{pseudoautosomal.hg19}}, and
   \code{\link[GWASTools]{pseudoautosomal.hg38}} from the
   \pkg{GWASTools} package (see details below).}
  \item{ignoreGaps}{skip assembly gaps only if larger than this
    threshold; in turn, if two unmasked regions are separated by an
    assembly gap not larger than \code{ignoreGaps}, they are joined in
    the resuling \code{\linkS4class{GRanges}} object.}
  \item{activeMasks}{masks to apply for determining unmasked region;
    defaults to the masks that are
    active by default in the \code{\link[BSgenome:BSgenome]{MaskedBSgenome}}
    object \code{x}. Therefore, this argument only needs to be set if a
    masking other than the default is necessary.}
}
\details{This function takes a
  \code{\link[BSgenome:BSgenome]{MaskedBSgenome}} object \code{x} and
  extracts the genomic regions that are unmasked in this genome, where
  the set of masks to apply can be specified using the
  \code{activeMasks} argument. The result is returned as a
  \code{\linkS4class{GRangesList}} object each component of which
  corresponds to one chromosome of the genome \code{x} - or a subset
  thereof if the \code{chrs} argument has been specified.

  The \code{pseudoautosomal} argument allows for a special treatment of
  pseudoautosomal regions. If not \code{NULL}, this argument must be
  a data frame that contains columns with names \dQuote{chrom}, 
  \dQuote{start.base}, and \dQuote{end.base}. The \dQuote{chrom}
  column must contain chromosome names as they appear in the
  \code{\link[BSgenome:BSgenome]{MaskedBSgenome}} object \code{x}. 
  The columns \dQuote{start.base} and \dQuote{end.base} must
  contain numeric values that specify the starts and ends of
  pseudoautosomal regions, respectively. The function is implemented
  such that the data frames \code{\link[GWASTools]{pseudoautosomal.hg18}},
  \code{\link[GWASTools]{pseudoautosomal.hg19}}, and
  \code{\link[GWASTools]{pseudoautosomal.hg38}} provided by the
  \pkg{GWASTools}
  package can be used (except for the chromosome names that need to be
  adapted to hg18/hg19/hg38). If the \code{pseudoautosomal} argument is
  specified correctly, the \code{unmaskedRegions} function produces
  separate components in the resulting \code{\linkS4class{GRangesList}}
  object - one for each pseudoautosomal region. These components are
  named as the corresponding row names in the data frame
  \code{pseudoautosomal}. Moreover, these regions are omitted from the
  list of unmasked regions of the chromosomes they are on.
}
\value{
  a \code{\linkS4class{GRangesList}} object (see details above)
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\linkS4class{GRangesList}},
  \code{\link[GWASTools:pseudoautosomal]{pseudoautosomal}}}
\examples{
## load packages to obtain masked hg38genome and
##  pseudoautosomal.hg19 from GWASTools package
if (require(BSgenome.Hsapiens.UCSC.hg38.masked) && require(GWASTools))
{
    ## extract unmasked regions of all autosomal chromosomes
    regions <- unmaskedRegions(BSgenome.Hsapiens.UCSC.hg38.masked,
                               chrs=paste0("chr", 1:22))
    names(regions)
    regions$chr1

    ## adjust chromosome names
    pseudoautosomal.hg38
    psaut <- pseudoautosomal.hg38
    psaut$chrom <- paste0("chr", psaut$chrom)
    psaut

    ## extract unmasked regions of sex chromosomes taking pseudoautosomal
    ## regions into account
    regions <- unmaskedRegions(BSgenome.Hsapiens.UCSC.hg38.masked,
                               chrs=c("chrX", "chrY"), pseudoautosomal=psaut)
    names(regions)
    regions$chrX
    regions$X.PAR1

    ## check overlap between X chromosome and a pseudoautosomal region
    intersect(regions$chrX, regions$X.PAR1)
}}
