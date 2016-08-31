\name{readRegionsFromBedFile}
\alias{readRegionsFromBedFile}
\title{Read Genomic Regions from BED File}
\description{
  Reads a BED file and returns the genomic regions as \code{GRanges} object
}
\usage{
readRegionsFromBedFile(file, header=FALSE, sep="\t",
                       col.names=c("chrom", "chromStart",
                                   "chromEnd", "names", "width",
                                   "strand"),
                       ignoreMcols=TRUE, seqInfo=NULL)
}
\arguments{
  \item{file}{the name of the file, text-mode connection, or URL to read
    data from}
  \item{header,sep,col.names}{arguments passed on to \code{\link{read.table}}}
  \item{ignoreMcols}{if \code{TRUE} (default), further columns are
    ignored; if \code{FALSE}, further columns are appended to the
    resulting \code{\linkS4class{GRanges}} object as metadata colums
    (see details below).}
  \item{seqInfo}{can be \code{NULL} (default) or an object of class
    \code{\linkS4class{Seqinfo}} (see details below).}
}
\details{This function is a simple wrapper around the
  \code{\link{read.table}} function that reads from a BED file
  and returns the genomic regions as a \code{\linkS4class{GRanges}}
  object. How the file is split into columns can be controlled by the
  arguments \code{header}, \code{sep}, and \code{col.names}. These
  arguments are passed on to \code{\link{read.table}} as they are.
  The choice of the \code{col.names} argument is crucial. A wrong
  \code{col.names} argument results in erroneous assignment of columns.
  The function \code{readRegionsFromBedFile} requires columns named
  \dQuote{chrom}, \dQuote{chromStart}, and \dQuote{chromEnd} to be
  present in the object returned from \code{\link{read.table}} upon
  reading from the BED file. If a column named \dQuote{strands} is
  contained in the BED file, this column is used as strand info in
  the resulting \code{\linkS4class{GRanges}} object.

  If \code{ignoreMcols=TRUE} (default), further columns are
  ignored. If \code{ignoreMcols=FALSE}, all columns other than
  \dQuote{chrom}, \dQuote{chromStart}, \dQuote{chromEnd},
  \dQuote{names}, \dQuote{strand}, and \dQuote{width} are appended
  to the resulting \code{\linkS4class{GRanges}} object as metadata
  columns.

  The \code{seqInfo} argument can be used to assign the right metadata,
  such as, genome, chromosome names, and chromosome lengths
  to the resulting \code{\linkS4class{GRanges}} object.
}
\value{
  a \code{\linkS4class{GRanges}} object
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}

\url{http://genome.ucsc.edu/FAQ/FAQformat.html#format1}
}
\seealso{\code{\link{read.table}}}
\examples{
## basic example (hg38 regions of HBA1 and HBA2)
bedFile <- system.file("examples/HBA.bed", package="podkat")
readRegionsFromBedFile(bedFile)

## example with enforcing seqinfo
data(hg38Unmasked)
readRegionsFromBedFile(bedFile, seqInfo=seqinfo(hg38Unmasked))

##
## example with regions targeted by Illumina TruSeq Exome Enrichment kit:
## download file "truseq_exome_targeted_regions.hg19.bed.chr.gz" from
## http://support.illumina.com/downloads/truseq_exome_targeted_regions_bed_file.ilmn
## (follow link "TruSeq Exome Targeted Regions BED file"; these regions
##  are based on hg19)
##
\dontrun{
readRegionsFromBedFile("truseq_exome_targeted_regions.hg19.bed.chr.gz")

data(hg19Unmasked)
readRegionsFromBedFile("truseq_exome_targeted_regions.hg19.bed.chr.gz",
                       seqInfo=seqinfo(hg19Unmasked))}
}
