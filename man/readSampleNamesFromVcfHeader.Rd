\name{readSampleNamesFromVcfHeader}
\alias{readSampleNamesFromVcfHeader}
\title{Read Sample Names from VCF File Header}
\description{
  Reads the header of a VCF file and returns sample names as character vector
}
\usage{
readSampleNamesFromVcfHeader(file, ...)
}
\arguments{
  \item{file}{a \code{\linkS4class{TabixFile}} object or a character
    string with a file name of the VCF file to read from; if \code{file}
    is a file name, the method internally creates a
    \code{\linkS4class{TabixFile}} object for this file name.}
  \item{...}{all additional arguments are passed on internally to
    \code{\link{scanBcfHeader}} function from the
    \code{\link[Rsamtools:Rsamtools-package]{Rsamtools}} package.}
}
\details{This function is a simple wrapper around the
  \code{\link{scanBcfHeader}} function from the
  \code{\link[Rsamtools:Rsamtools-package]{Rsamtools}} package
  that scans the header of a VCF file and returns the sample
  names as a character vector.
}
\value{
  a character vector with sample names
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}

\url{http://www.1000genomes.org/wiki/analysis/variant-call-format/vcf-variant-call-format-version-42}

Li, H., Handsaker, B., Wysoker, A., Fenell, T., Ruan, J., Homer, N.,
Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data
Processing Subgroup (2009) The Sequence Alignment/Map format and
SAMtools. \emph{Bioinformatics} \bold{25}, 2078-2079.
}
\seealso{\code{\link{scanBcfHeader}}}
\examples{
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
readSampleNamesFromVcfHeader(vcfFile)
}
