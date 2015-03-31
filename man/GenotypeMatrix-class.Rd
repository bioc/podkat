\name{GenotypeMatrix-class}
\docType{class}
\alias{GenotypeMatrix-class}
\alias{class:GenotypeMatrix}
\alias{GenotypeMatrix}
\alias{[,GenotypeMatrix,index,missing,missing-method}
\alias{[,GenotypeMatrix,missing,index,missing-method}
\alias{[,GenotypeMatrix,index,index,missing-method}
\alias{variantInfo,GenotypeMatrix-method}
\alias{MAF,GenotypeMatrix-method}
\alias{show,GenotypeMatrix-method}

\title{Class \code{GenotypeMatrix}}
\description{S4 class for storing genotypes efficiently as
  column-oriented sparse matrices along with variant info}
\section{Details}{
  This class stores genotypes as a column-oriented sparse numeric
  matrix, where rows correspond to samples and columns correspond to
  variants. This is accomplished by extending the
  \code{\linkS4class{dgCMatrix}} class from which this class inherits
  all slots. Information about variants is stored in an additional slot
  named \code{variantInfo}. This slot must be of class
  \code{\linkS4class{VariantInfo}} and have exactly as many elements as the
  genotype matrix has columns. The \code{variantInfo} slot has a
  dedicated metadata column named \dQuote{MAF} that contains the minor
  allele frequencies (MAFs) of the variants. For convenience,
  accessor functions \code{variantInfo} and \code{MAF} are available (see
  below).
  
  Objects of this class should only be created and manipulated by the
  constructors and accessors described below, as only these methods ensure the
  integrity of the created objects. Direct modification of object slots
  is strongly discouraged!
}
\section{Constructors}{
  See help pages \code{\link{genotypeMatrix}} and
  \code{\link{readGenotypeMatrix}}.
}
\section{Methods}{
   \describe{
     \item{show}{\code{signature(object="GenotypeMatrix")}:
       displays the matrix dimensions (i.e. the number of samples and variants)
       along with some basic statistics of the minor allele frequency
       (MAF).}
     }
}
\section{Accessors}{
   \describe{
     \item{variantInfo}{\code{signature(object="GenotypeMatrix")}:
       returns variant information as a \code{\linkS4class{VariantInfo}}
       object.}
     \item{MAF}{\code{signature(object="GenotypeMatrix")}:
       returns a numeric vector with the minor allele frequencies (MAFs).}
   }
   Row and column names can be set and get as usual for matrix-like
   objects with \code{\link[base:rownames]{rownames}} and
   \code{\link[base:colnames]{colnames}}, respectively.
   When setting the column names of a \code{GenotypeMatrix}
   object, both the names of the variant info (slot \code{variantInfo}) and the
   column names of the matrix are set.
}
\section{Subsetting}{
  In the following code snippets, \code{x} is a
  \code{\linkS4class{GenotypeMatrix}} object.
  \describe{
    \item{}{\code{x[i, ]}: returns a
      \code{\linkS4class{GenotypeMatrix}} object that only
      contains the samples selected by the index vector \code{i}
    }
    \item{}{\code{x[, j]}: returns a
      \code{\linkS4class{GenotypeMatrix}}
      object that only contains the variants selected by the index
      vector \code{j}
    }
    \item{}{\code{x[i, j]}: returns a
      \code{\linkS4class{GenotypeMatrix}}
      object that only contains the samples selected by the index vector
      \code{i} and the variants selected by the index vector \code{j}
    }
  }
  None of these subsetting functions support a \code{drop} argument.
  As soon as a \code{drop} argument is supplied, no matter whether
  \code{TRUE} or \code{FALSE}, all variant information is stripped off and a
  \code{\linkS4class{dgCMatrix}} object is returned.

  By default, MAFs
  are not altered by subsetting samples. However, if the optional
  argument \code{recomputeMAF} is set to \code{TRUE} (the default is
  \code{FALSE}), MAFs are recomputed for the resulting
  subsetted genotype matrix as described in
  \code{\link{genotypeMatrix}}. The ploidy for computing MAFs can
  be controlled by the optional \code{ploidy} argument (the default is 2).
}

\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\linkS4class{dgCMatrix}}, \code{\linkS4class{VariantInfo}},
\code{\link{genotypeMatrix}}, \code{\link{readGenotypeMatrix}}}
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
\keyword{classes}
