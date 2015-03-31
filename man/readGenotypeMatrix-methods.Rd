\name{readGenotypeMatrix-methods}
\docType{methods}
\alias{readGenotypeMatrix-methods}
\alias{method:readGenotypeMatrix}
\alias{readGenotypeMatrix}
\alias{readGenotypeMatrix,TabixFile,GRanges-method}
\alias{readGenotypeMatrix,TabixFile,missing-method}
\alias{readGenotypeMatrix,character,GRanges-method}
\alias{readGenotypeMatrix,character,missing-method}
\title{Read from VCF File}
\description{A fast lightweight function that reads from a VCF file and
  returns the result as a \code{\linkS4class{GenotypeMatrix}} object}
\usage{
\S4method{readGenotypeMatrix}{TabixFile,GRanges}(file, regions, subset,
                   noIndels=TRUE, onlyPass=TRUE,
                   na.limit=1, MAF.limit=1,
                   na.action=c("impute.major", "omit", "fail"),
                   MAF.action=c("invert", "omit", "ignore", "fail"),
                   sex=NULL)
\S4method{readGenotypeMatrix}{TabixFile,missing}(file, regions, ...)
\S4method{readGenotypeMatrix}{character,GRanges}(file, regions, ...)
\S4method{readGenotypeMatrix}{character,missing}(file, regions, ...)
}
\arguments{
  \item{file}{a \code{\linkS4class{TabixFile}} object or a character
    string with a file name of the VCF file to read from; if \code{file}
    is a file name, the method internally creates a
    \code{\linkS4class{TabixFile}} object for this file name.}
  \item{regions}{a \code{\linkS4class{GRanges}} object that specifies
    which genomic regions to read from the VCF file; if missing,
    the entire VCF file is read.}
  \item{subset}{a numeric vector with indices or a character vector with
    names of samples to restrict to; if specified, only these samples'
    genotypes are read from the VCF file and all other samples are ignored
    and omitted from the \code{\linkS4class{GenotypeMatrix}} object that
    is returned. Moreover, minor allele frequencies (MAFs) are only
    computed from the genotypes of the samples specified by \code{subset}.}
  \item{noIndels}{if \code{TRUE} (default), only single-nucleotide
    variants (SNVs) are considered and indels are skipped.}
  \item{onlyPass}{if \code{TRUE} (default), only variants are considered
    whose value in the \code{FILTER} column is \dQuote{PASS}.}
   \item{na.limit}{all variants with a missing value ratio above this
    threshold will be omitted from the output object.}
  \item{MAF.limit}{all variants with an MAF above this
    threshold will be omitted from the output object.}
  \item{na.action}{if \dQuote{impute.major}, all missing values will
    be imputed by major alleles in the output object. If
    \dQuote{omit}, all variants containing missing values will be
    omitted in the output object. If \dQuote{fail}, the function stops
    with an error if a variant contains any missing values.}
  \item{MAF.action}{if \dQuote{invert}, all variants with an MAF
    exceeding 0.5 will be inverted in the sense that all minor alleles
    will be replaced by major alleles and vice versa. If
    \dQuote{omit}, all variants with an MAF greater
    than 0.5 are omitted in the output object. If \dQuote{ignore},
    no action is taken and MAFs greater than 0.5 are kept as they are.
    If \dQuote{fail}, the function stops with an error if any
    variant has an MAF greater than 0.5.}
  \item{sex}{if \code{NULL}, all samples are treated the same
    without any modifications; if \code{sex} is a factor with levels
    \code{F} (female) and \code{M} (male) that is as long as
    \code{subset} or as the VCF file has samples, this argument is
    interpreted as the sex of the samples. In this case, the genotypes
    corresponding to male samples are
    doubled before further processing. This is designed for mixed-sex
    analyses of the X chromosome outside of the pseudoautosomal
    regions.}
  \item{...}{for the three latter methods above, all other parameters
    are passed on to the method with signature \code{TabixFile,GRanges}.}
}    
\details{This method uses the \code{tabix} API provided by the
  \pkg{Rsamtools} package
  to read from a VCF file, parses the result into a sparse matrix
  along with positional information, and returns the result as a
  \code{\linkS4class{GenotypeMatrix}} object. Reading can be restricted
  to certain regions by specifying the \code{regions} object.
  Note that it might not be possible to read a very large VCF file
  as a whole. 

  For all variants, filters in terms of missing values and MAFs can be
  applied. Moreover, variants with MAFs greater than 0.5 can filtered
  out or inverted. For details, see descriptions of parameters
  \code{na.limit}, \code{MAF.limit}, \code{na.action}, and
  \code{MAF.action} above.
}
\value{returns an object of class \code{\linkS4class{GenotypeMatrix}}
}

\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}

\url{http://www.1000genomes.org/wiki/analysis/variant-call-format/vcf-variant-call-format-version-42}

Li, H., Handsaker, B., Wysoker, A., Fenell, T., Ruan, J., Homer, N.,
Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data
Processing Subgroup (2009) The Sequence Alignment/Map format and
SAMtools. \emph{Bioinformatics} \bold{25}, 2078-2079.
}
\seealso{\code{\linkS4class{GenotypeMatrix}}}
\examples{
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
readGenotypeMatrix(vcfFile)
readGenotypeMatrix(vcfFile, onlyPass=FALSE, MAF.action="ignore")
}
\keyword{methods}
