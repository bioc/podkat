\name{readVariantInfo-methods}
\docType{methods}
\alias{readVariantInfo-methods}
\alias{method:readVariantInfo}
\alias{readVariantInfo}
\alias{readVariantInfo,TabixFile,GRanges-method}
\alias{readVariantInfo,TabixFile,missing-method}
\alias{readVariantInfo,character,GRanges-method}
\alias{readVariantInfo,character,missing-method}
\title{Read information about variants from VCF file}
\description{A fast lightweight function that determines information on
  variants ocurring in a VCF file and
  returns the result as a \code{\linkS4class{VariantInfo}} object}
\usage{
\S4method{readVariantInfo}{TabixFile,GRanges}(file, regions, subset,
                noIndels=TRUE, onlyPass=TRUE,
                na.limit=1, MAF.limit=1,
                na.action=c("impute.major", "omit", "fail"),
                MAF.action=c("ignore", "omit","invert", "fail"),
                omitZeroMAF=TRUE, refAlt=FALSE, sex=NULL)
\S4method{readVariantInfo}{TabixFile,missing}(file, regions, ...)
\S4method{readVariantInfo}{character,GRanges}(file, regions, ...)
\S4method{readVariantInfo}{character,missing}(file, regions, ...)
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
    genotypes are considered when determining the minor allele
    frequencies (MAFs) of variants.}
  \item{noIndels}{if \code{TRUE} (default), only single-nucleotide
    variants (SNVs) are considered and indels are skipped.}
  \item{onlyPass}{if \code{TRUE} (default), only variants are considered
    whose value in the \code{FILTER} column is \dQuote{PASS}.}
   \item{na.limit}{all variants with a missing value ratio above this
    threshold will be omitted from the output object.}
  \item{MAF.limit}{all variants with an MAF above this
    threshold will be omitted from the output object.}
  \item{na.action}{if \dQuote{impute.major}, all missing values are
    considered as major alleles when computing MAFs. If
    \dQuote{omit}, all variants containing missing values will be
    omitted in the output object. If \dQuote{fail}, the function stops
    with an error if a variant contains any missing values.}
  \item{MAF.action}{if \dQuote{ignore} (default),
    no action is taken for variants with an MAF greater than 0.5, these
    variants are kept and included in the output object as they are.
    If \dQuote{omit}, all variants with an MAF greater
    than 0.5 are omitted in the output object. If \dQuote{fail}, the
    function stops with an error if any variant has an MAF greater than
    0.5. If \dQuote{invert}, all variants with an MAF
    exceeding 0.5 will be inverted in the sense that all minor alleles
    will be replaced by major alleles and vice versa. Note: if this setting is
    used in conjunction with \code{refAlt=TRUE}, the MAFs of the
    variants that have been inverted do no longer correspond to the
    true alternate allele.}
  \item{omitZeroMAF}{if \code{TRUE} (default), variants with an MAF
    of 0 are not considered and omitted from the output object.}
  \item{refAlt}{if \code{TRUE}, two metadata columns named \dQuote{ref} and
    \dQuote{alt} are added to the output object that contain
    reference and alternate alleles. Note that these sequences can be
    quite long for indels, which may result in large memory consumption.
    The default is \code{FALSE}.}
  \item{sex}{if \code{NULL}, all samples are treated the same
    without any modifications; if \code{sex} is a factor with levels
    \code{F} (female) and \code{M} (male) that is as long as
    \code{subset} or as the VCF file has samples, this argument is
    interpreted as the sex of the samples. In this case, the genotypes
    corresponding to male samples are doubled before computing MAFs.
    The option to supply the \code{sex} argument is meant to allow
    for a correct estimate of MAFs as \code{\link{readGenotypeMatrix}}
    and \code{\link{assocTest}} compute it. Note, however, that the
    MAFs computed in this way do not correspond to the true MAFs
    contained in the data.}
  \item{...}{for the three latter methods above, all other parameters
    are passed on to the method with signature
    \code{TabixFile,GRanges}.}
}    
\details{This method uses the \dQuote{tabix} API provided by the
  \pkg{Rsamtools} package to parse a VCF file.
  The \code{readVariantInfo} method considers each variant and
  determines its minor allele frequency (MAF) and the type of the
  variant. The result is returned as a \code{\linkS4class{VariantInfo}}
  object, i.e. a \code{\linkS4class{GRanges}} object with two
  metadata columns \dQuote{MAF} and \dQuote{type}. The former contains
  the MAF of each variant, while the latter is a factor column that
  contains information about the type of the variant. Possible values
  in this column are \dQuote{INDEL} (insertion or deletion),
  \dQuote{MULTIPLE} (single-nucleotide variant with multiple alternate
  alleles), \dQuote{TRANSITION} (single-nucleotide variation A/G or C/T),
  \dQuote{TRANSVERSION} (single-nucleotide variation A/C, A/T, C/G, or
  G/T), or \dQuote{UNKNOWN} (anything else). If \code{refAlt} is \code{TRUE},
  two further metadata columns \dQuote{ref} and \dQuote{alt} are
  included which contain reference and alternate alleles of each variant.

  For all variants, filters in terms of missing values and MAFs can be
  applied. Moreover, variants with MAFs greater than 0.5 can filtered
  out or inverted. For details, see descriptions of parameters
  \code{na.limit}, \code{MAF.limit}, \code{na.action}, and
  \code{MAF.action} above.
}
\value{returns an object of class \code{\linkS4class{VariantInfo}}
}

\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}

\url{https://github.com/samtools/hts-specs}

Li, H., Handsaker, B., Wysoker, A., Fenell, T., Ruan, J., Homer, N.,
Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data
Processing Subgroup (2009) The Sequence Alignment/Map format and
SAMtools. \emph{Bioinformatics} \bold{25}, 2078-2079.
DOI: \doi{10.1093/bioinformatics/btp352}.
}
\seealso{\code{\linkS4class{GenotypeMatrix}}}
\examples{
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")

## default parameters
vInfo <- readVariantInfo(vcfFile)
vInfo
summary(vInfo)

## including zero MAF variants and reference/alternate alleles
vInfo <- readVariantInfo(vcfFile, omitZeroMAF=FALSE, refAlt=TRUE)
vInfo
summary(vInfo)
}
\keyword{methods}
