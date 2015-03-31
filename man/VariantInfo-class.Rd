\name{VariantInfo-class}
\docType{class}
\alias{VariantInfo-class}
\alias{class:VariantInfo}
\alias{VariantInfo}
\alias{variantInfo}
\alias{variantInfo,missing-method}
\alias{variantInfo,GRanges-method}
\alias{MAF}
\alias{MAF,VariantInfo-method}
\alias{show,VariantInfo-method}
\alias{summary,VariantInfo-method}

\title{Class \code{VariantInfo}}
\description{S4 class for storing information about variants}
\section{Details}{
  This class extends the class \code{\linkS4class{GRanges}} without
  adding any extra slots. The main difference is that \code{VariantInfo}
  objects always have a metadata column \dQuote{MAF} that contains
  minor allele frequencies (MAFs). A special
  \code{summary} method allows for computing statistics about
  MAFs and types of variants.
 
  Objects of this class should only be created and manipulated by the
  constructors and accessors described below, as only these methods ensure the
  integrity of the created objects. Direct modification of object slots
  is strongly discouraged!
}
\section{Constructors}{
  \describe{
    \item{variantInfo}{\code{signature(x="missing")}:
      creates an empty \code{VariantInfo} object}
    \item{variantInfo}{\code{signature(x="GRanges")}:
      coerces a \code{\linkS4class{GRanges}} object to class
      \code{VariantInfo} by adding a \dQuote{MAF} metadata column that
      is initialized with \code{NA} values.}
  }
  Furthermore, see the help page of \code{\link{readVariantInfo}}.
}
\section{Accessors}{
   \describe{
     \item{MAF}{\code{signature(object="VariantInfo")}:
       returns a numeric vector with the minor allele frequencies (MAFs).}
   }
}
\section{Methods}{
   \describe{
     \item{summary}{\code{signature(object="VariantInfo")}:
       computes and prints a summary about the MAFs and variant types
       (if available); this method returns a list with summarized values.}
   }
   All other methods, including sub-setting, are inherited from the
   \code{\linkS4class{GRanges}} class.
}

\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\linkS4class{GRanges}}, \code{\link{readVariantInfo}}, 
\code{\link{genotypeMatrix}}, \code{\link{readGenotypeMatrix}}}
\examples{
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")

## read variant info directly from VCF file
vInfo <- readVariantInfo(vcfFile, omitZeroMAF=FALSE)
vInfo
summary(vInfo)

## read full genotype from VCF file
geno <- readGenotypeMatrix(vcfFile)

## display summary of variant info stored in genotype matrix
summary(variantInfo(geno))
}
\keyword{classes}
