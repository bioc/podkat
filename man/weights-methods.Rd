\name{weights}
\docType{methods}
\alias{weights}
\alias{weights-methods}
\alias{weights,AssocTestResult-method}
\alias{weights,AssocTestResultRanges-method}
\title{Extract Contribution Weights of Variants}
\description{
  Method for extracting the contributions that each variant makes to the
  test statistic of an association test
}
\usage{
\S4method{weights}{AssocTestResult}(object, Z, model)
\S4method{weights}{AssocTestResultRanges}(object, Z, model, limit=20, sex=NULL)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{AssocTestResult}}
    or \code{\linkS4class{AssocTestResultRanges}}}
  \item{Z}{an object of class \code{\linkS4class{GenotypeMatrix}},
    an object of class \code{\link[Rsamtools]{TabixFile}}, or a
    character string with the file name of a VCF file}
  \item{model}{an object of class \code{\linkS4class{NullModel}}}
  \item{limit}{maximum number of regions to be processed; set to
    \code{Inf} or non-numeric value like \code{NA} or \code{NULL}
    to disable limitation. Do this with caution, in particular, when
    reading from a VCF file, as reading of excessively large regions
    from VCF files may take very long or even kill the R session because
    of excessive memory comsumption!}
  \item{sex}{if \code{NULL}, all samples are treated the same
    without any modifications; if \code{sex} is a factor with levels
    \code{F} (female) and \code{M} (male) that is as long as
    the number of samples in \code{model}, this argument is
    interpreted as the sex of the samples. In this case, the genotypes
    corresponding to male samples are
    doubled before further processing. This is designed for mixed-sex
    analyses of the X chromosome outside of the pseudoautosomal
    regions.}
}
\details{
  Upon successful completion of an association test, the \code{weights}
  method allows for finding out the individual contributions each of the
  variants made to the test statistic. This computation is only possible
  for kernels \dQuote{linear.podkat} and \dQuote{linear.SKAT}
  (see \code{\link{computeKernel}}).

  If called for an \code{\linkS4class{AssocTestResult}} object as
  first argument \code{object}, a \code{\linkS4class{GenotypeMatrix}}
  object \code{Z}, and a \code{\linkS4class{NullModel}} object
  \code{model}, \code{weights} returns a
  \code{\linkS4class{GRanges}} object that contains all variants of
  \code{variantInfo(Z)} along with two numerical metadata columns named
  \dQuote{weight.raw} and \dQuote{weight.contribution}. The column
  \dQuote{weight.raw} corresponds to raw contributions. These are
  signed, i.e. a positive value indicates a positive association, while
  a negative value indicates a negative association. The larger the
  absolute value, the larger the contribution. The column
  \dQuote{weight.contribution} corresponds to relative contributions.
  These values are non-negative and they sum up to 1.
  For mathematical details, see Subsection 9.4 of the package vignette.

  If \code{weights} is called for an
  \code{\linkS4class{AssocTestResultRanges}} object \code{object},
  a second argument \code{Z} that is an object of class
  \code{\linkS4class{GenotypeMatrix}}, an object of class
  \code{\link[Rsamtools]{TabixFile}}, or a character string with the
  name of a VCF file, and a \code{\linkS4class{NullModel}} object
  \code{model}, the contribution weights described above
  are computed for each region in \code{object}. In this case, the
  method returns a \code{\linkS4class{GRangesList}} with as many
  components as \code{object} has regions, where each list component is
  a \code{\linkS4class{GRanges}} object containing the
  contribution weights as described above.

  It is essential for \code{weights} to work correctly that
  \code{object} is actually the result of an association test between
  \code{Z} and \code{model}. If called for objects that actually do not
  belong to each other, the results are void. The method is implemented
  such that all possible checks are made that possibly detect
  inconsistencies between the input objects. However, the final
  responsibility is left to user to make sure that all data are
  consistent. Special caution is necessary if \code{weights} is run
  for an \code{\link{AssocTestResultRanges}} object that has
  been obtained by merging multiple \code{\linkS4class{AssocTestResultRanges}}
  using the \code{c} method. The \code{c} method performs several checks
  to ensure consistency of association test parameters among the merged
  results, but the \code{sex} parameter is an exception: if it appears
  to be inconsistent among the results to merge, it is omitted from the
  merged object (see also \code{\linkS4class{AssocTestResultRanges}}).

  The \code{weights} method needs to re-evaluate some computations of
  the association test. In case it is called for \code{Z} being a
  \code{\link[Rsamtools]{TabixFile}} object or file name of a VCF file,
  \code{weights} even needs to re-read the genotype data from the file.
  Therefore, the method has a safety limit not to process too many
  regions (see \code{limit} argument described above).
}
\value{
  an object of class \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} (see details above)
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}
}
\seealso{\code{\link{assocTest}}, \code{\link{AssocTestResult}},
  \code{\link{AssocTestResultRanges}}, \code{\link{nullModel}},
  \code{\linkS4class{NullModel}}, \code{\link{computeKernel}},
  \code{\linkS4class{GenotypeMatrix}}, \code{\link{p.adjust}},
  \code{\link{filterResult}}}
\examples{
## load genome description
data(hgA)

## partition genome into overlapping windows
windows <- partitionRegions(hgA)

## load genotype data from VCF file
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
Z <- readGenotypeMatrix(vcfFile)

## read phenotype data from CSV file (continuous trait + covariates)
phenoFile <- system.file("examples/example1lin.csv", package="podkat")
pheno <- read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
model <- nullModel(y ~ ., pheno)

## perform association test
res <- assocTest(Z, model, windows)

## perform multiple testing correction and filter for
## significant regions
res <- filterResult(p.adjust(res), filterBy="p.value.adj")

## compute contributions
contrib <- weights(res, Z, model)
contrib

## extract most indicative variants
filterResult(contrib)

## plot contributions
plot(contrib[[1]], "weight.raw")
plot(contrib[[1]], "weight.contribution", type="b", alongGenome=TRUE)
}
\keyword{methods}
