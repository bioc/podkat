\name{genotypeMatrix-methods}
\docType{methods}
\alias{genotypeMatrix-methods}
\alias{method:genotypeMatrix}
\alias{genotypeMatrix}
\alias{genotypeMatrix,ANY,GRanges,missing-method}
\alias{genotypeMatrix,ANY,numeric,character-method}
\alias{genotypeMatrix,ANY,character,missing-method}
\alias{genotypeMatrix,ANY,missing,missing-method}
\alias{genotypeMatrix,eSet,numeric,character-method}
\alias{genotypeMatrix,eSet,character,missing-method}
\alias{genotypeMatrix,eSet,character,character-method}

\title{Constructors for Creating \code{\linkS4class{GenotypeMatrix}} Objects}
\description{Create \code{\linkS4class{GenotypeMatrix}} object from
  (sparse) matrix object and positions of variants}
\usage{\S4method{genotypeMatrix}{ANY,GRanges,missing}(Z, pos, seqnames,
       ploidy=2, na.string=NULL, na.limit=1, MAF.limit=1,
       na.action=c("impute.major", "omit", "fail"),
       MAF.action=c("invert", "omit", "ignore", "fail"),
       sex=NULL)
\S4method{genotypeMatrix}{ANY,numeric,character}(Z, pos, seqnames, ...)
\S4method{genotypeMatrix}{ANY,character,missing}(Z, pos, seqnames, ...)
\S4method{genotypeMatrix}{ANY,missing,missing}(Z, pos, seqnames, subset,
       noIndels=TRUE, onlyPass=TRUE, sex=NULL, ...)
\S4method{genotypeMatrix}{eSet,numeric,character}(Z, pos, seqnames, ...)
\S4method{genotypeMatrix}{eSet,character,missing}(Z, pos, seqnames, ...)
\S4method{genotypeMatrix}{eSet,character,character}(Z, pos, seqnames, ...)
}
\arguments{
  \item{Z}{an object of class \code{\linkS4class{dgCMatrix}}, a numeric
    matrix, a character matrix, an object of class
    \code{VCF}, or an object of class \code{\link[Biobase]{eSet}}
    (see details below)}
  \item{pos}{an object of class \code{\linkS4class{GRanges}}, a numeric
    vector, or a character vector (see details below)}
  \item{seqnames}{a character vector (see details below)}
  \item{ploidy}{determines the ploidy of the genome for the computation of
    minor allele frequencies (MAFs) and the possible inversion of
    columns with an MAF exceeding
    0.5; the elements of \code{Z} may not exceed this value.}
  \item{subset}{a numeric vector with indices or a character vector with
    names of samples to restrict to}
  \item{na.limit}{all columns with a missing value ratio above this
    threshold will be omitted from the output object.}
  \item{MAF.limit}{all columns with an MAF above this
    threshold will be omitted from the output object.}
  \item{na.action}{if \dQuote{impute.major}, all missing values will
    be imputed by major alleles in the output object. If
    \dQuote{omit}, all columns containing missing values will be
    omitted in the output object. If \dQuote{fail}, the function stops
    with an error if \code{Z} contains any missing values.}
  \item{MAF.action}{if \dQuote{invert}, all columns with an MAF
    exceeding 0.5 will be inverted in the sense that all minor alleles
    will be replaced by major alleles and vice versa. For numerical
    \code{Z}, this is accomplished by subtracting the column from the
    ploidy value. If \dQuote{omit}, all columns with an MAF greater
    than 0.5 are omitted in the output object. If \dQuote{ignore},
    no action is taken and MAFs greater than 0.5 are kept as they are.
    If \dQuote{fail}, the function stops with an error if \code{Z}
    contains any column with an MAF greater than 0.5.}
  \item{noIndels}{if \code{TRUE} (default), only single nucleotide
    variants (SNVs) are considered
    and indels are skipped; only works if the \code{ALT} column is
    present in the \code{VCF} object \code{Z}, otherwise a
    warning is shown and the \code{noIndels} argument is ignored.}
  \item{onlyPass}{if \code{TRUE} (default), only variants are considered
    whose value in the \code{FILTER} column is \dQuote{PASS};
    only works if the \code{FILTER} column is
    present in the \code{VCF} object \code{Z}, otherwise a
    warning is shown and the \code{onlyPass} argument is ignored.}
  \item{na.string}{if not \code{NULL}, all \dQuote{.} entries in the character
    matrix or \code{VCF} genotype are replaced with this
    string before parsing the matrix.}
  \item{sex}{if \code{NULL}, all rows of \code{Z} are treated the same
    without any modifications; if \code{sex} is a factor with levels
    \code{F} (female) and \code{M} (male) that is as long as
    \code{Z} has rows, this argument is interpreted as the sex of the
    samples. In this case, the rows corresponding to male samples are
    doubled before further processing. This is designed for mixed-sex
    analyses of the X chromosome outside of the pseudoautosomal
    regions.}
  \item{...}{all additional arguments are passed on internally to the
    \code{genotypeMatrix} method with signature \code{ANY,GRanges,missing}.}
}    
\details{This method provides different ways of constructing an object of
  class \code{\linkS4class{GenotypeMatrix}} from other types of objects.
  The typical case is when a matrix object is combined with positional
  information. The first three variants listed above work with \code{Z}
  being a \code{\linkS4class{dgCMatrix}} object, a numeric matrix, or
  a character matrix.

  If \code{Z} is a \code{\linkS4class{dgCMatrix}} object or a matrix,
  rows are interpreted as samples and columns are interpreted as variants.
  For \code{\linkS4class{dgCMatrix}} objects and numeric
  matrices, matrix entries are interpreted as the numbers of minor alleles
  (with 0 meaning only major alleles). In this case, minor allele
  frequencies (MAFs) are computed as column sums divided by the number
  of alleles, i.e. the number of samples/rows multiplied by the
  \code{ploidy} parameter. If \code{Z} is a character matrix, the matrix
  entries need to comply to the format of the \dQuote{GT} field in
  VCF files. MAFs are computed as the actual relative frequency of
  minor alleles among all alleles in a column. For a diploid genome,
  therefore, this results in the same MAF estimate as mentioned above.
  However, some VCF readers, most importantly \code{readVcf} from the
  \pkg{VariantAnnotation} package, replace missing genotypes by a
  single \dQuote{.} even for non-haploid genomes, which would result in a
  wrong MAF estimate. To correct for this, the \code{na.string}
  parameter is available. If not \code{NULL}, all \dQuote{.} entries
  in the matrix are replaced by \code{na.string} before parsing the
  matrix. The correct setting for a diploid genome would be
  \dQuote{./.}.

  Positional information can be passed to the function in three
  different ways:
  \itemize{
    \item{by supplying a \code{\linkS4class{GRanges}} object as
      \code{pos} argument and omitting the \code{seqnames} argument,}
    \item{by supplying a numeric vector of positions as \code{pos}
      argument and sequence/chromosome names as \code{seqnames}
      argument, or}
    \item{by supplying a character vector with entries of the format
      \dQuote{seqname:pos} as \code{pos} argument and omitting the
      \code{seqnames} argument.}
  }
  In all three cases, the lengths of the arguments \code{pos} and
  \code{seqnames} (if not omitted) must match the number of columns of
  \code{Z}.

  If the arguments \code{pos} and \code{seqnames} are not specified,
  argument \code{Z} can (and must) be an object of class
  \code{VCF} (cf. package \pkg{VariantAnnotation}). In this case,
  the \code{genotypeMatrix} method extracts both the genotype matrix and
  positional information directly from the \code{VCF}
  object. Consequently, the \code{VCF} object \code{Z} must contain genotype
  information. If so, the genotype matrix is parsed and converted as
  described above for character matrices. Moreover, indels and variants
  that did not pass all quality filters can be
  skipped (see description of arguments \code{noIndels} and
  \code{onlyPass} above).

  For all variants, filters in terms of missing values and MAFs can be
  applied. Moreover, variants with MAFs greater than 0.5 can filtered
  out or inverted. For details, see descriptions of parameters
  \code{na.limit}, \code{MAF.limit}, \code{na.action}, and
  \code{MAF.action} above.

  For convenience, \code{genotypeMatrix} also allows for converting
  SNP genotype matrices stored in \code{\link[Biobase]{eSet}} objects,
  e.g. \code{\link[Biobase]{SnpSet}} objects or \code{SnpSetIllumina}
  objects (cf. package \pkg{beadarraySNP}). If \code{genotypeMatrix} is
  called with an \code{\link[Biobase]{eSet}} object as first argument
  \code{Z}, the method first checks whether there is a slot \code{call}
  in \code{assayData(Z)} and whether it is a matrix. If so, this matrix
  is interpreted as follows: 1 corresponds to genotype \dQuote{AA},
  2 corresponds to the genotype \dQuote{Aa}, and 3 corresponds to the
  genotype \dQuote{aa}, where \dQuote{A} is the major allele and
  \dQuote{a} is the minor allele. If \code{pos} is a numeric vector and
  \code{seqnames} is a character vector or if \code{pos} is a character
  vector and \code{seqnames} is missing, then these two arguments
  are interpreted as described above. However, if \code{pos}
  and \code{seqnames} are both single strings (character vectors of
  length 1), then \code{pos} is interpreted as the name of the
  feature data column that contains positional information and
  \code{seqnames} is interpreted as the feature data column that
  contains the chromosome on which each variant is located.
  Correspondingly, \code{featureData(Z)[[pos]]} must be available
  and must be a numeric vector. Correspondingly,
  \code{featureData(Z)[[seqnames]]} must be available and must be
  a character vector (or a data type that can be cast to a character
  vector).
}
\value{returns an object of class \code{\linkS4class{GenotypeMatrix}}
}

\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}
  
\url{https://github.com/samtools/hts-specs}

Obenchain, V., Lawrence, M., Carey, V., Gogarten, S., Shannon, P., and
Morgan, M. (2014) VariantAnnotation: a Bioconductor package for
exploration and annotation of genetic variants. \emph{Bioinformatics}
\bold{30}, 2076-2078. DOI: \doi{10.1093/bioinformatics/btu168}.
}
\seealso{\code{\linkS4class{GenotypeMatrix}}, 
  \code{\linkS4class{dgCMatrix}}, \code{\linkS4class{GRanges}}
}
\examples{
## create a toy example
A <- matrix(rbinom(50, 2, prob=0.2), 5, 10)
sA <- as(A, "dgCMatrix")
pos <- sort(sample(1:10000, ncol(A)))
seqname <- "chr1"

## variant with 'GRanges' object
gr <- GRanges(seqnames=seqname, ranges=IRanges(start=pos, width=1))
gtm <- genotypeMatrix(A, gr)
gtm
as.matrix(gtm)
variantInfo(gtm)
MAF(gtm)

## variant with 'pos' and 'seqnames' object
genotypeMatrix(sA, pos, seqname)

## variant with 'seqname:pos' strings passed through 'pos' argument
spos <- paste(seqname, pos, sep=":")
spos
genotypeMatrix(sA, spos)

## read data from VCF file using 'readVcf()' from the 'VariantAnnotation'
## package
if (require(VariantAnnotation))
{
    vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
    sp <- ScanVcfParam(info=NA, genome="GT", fixed=c("ALT", "FILTER"))
    vcf <- readVcf(vcfFile, genome="hgA", param=sp)
    rowRanges(vcf)

    ## call constructor for 'VCF' object
    gtm <- genotypeMatrix(vcf)
    gtm
    variantInfo(gtm)

    ## alternatively, extract information from 'VCF' object and use
    ## variant with character matrix and 'GRanges' positions
    ## note that, in 'VCF' objects, rows correspond to variants and
    ## columns correspond to samples, therefore, we have to transpose the
    ## genotype
    gt <- t(geno(vcf)$GT)
    gt[1:5, 1:5]
    gr <- rowRanges(vcf)
    gtm <- genotypeMatrix(gt, gr)
    as.matrix(gtm[1:20, 1:5, recomputeMAF=TRUE])
}}
\keyword{methods}
