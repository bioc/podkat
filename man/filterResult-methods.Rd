\name{filterResult-methods}
\docType{methods}
\alias{filterResult}
\alias{filterResult-methods}
\alias{filterResult,AssocTestResultRanges-method}
\alias{filterResult,GRanges-method}
\alias{filterResult,GRangesList-method}
\title{Filter Association Test Results According to p-Values or
  Variants' Contributions}
\description{
  Given an \code{\linkS4class{AssocTestResultRanges}} object,
  this method filters regions according to p-values or variants' contributions.
}
\usage{
\S4method{filterResult}{AssocTestResultRanges}(object, cutoff=0.05,
             filterBy=c("p.value", "p.value.adj", "p.value.resampled",
             "p.value.resampled.adj"))
\S4method{filterResult}{GRanges}(object, cutoff=0.1)
\S4method{filterResult}{GRangesList}(object, cutoff=0.1)
}
\arguments{
  \item{object}{object of class
    \code{\linkS4class{AssocTestResultRanges}},
    \code{\linkS4class{GRanges}}, or \code{\linkS4class{GRangesList}}}
  \item{cutoff}{threshold}
  \item{filterBy}{according to which p-value column filtering should be
    done; the default is \dQuote{p.value}.}
}
\details{
  If called for an \code{\linkS4class{AssocTestResultRanges}} object as
  first argument, this
  method strips off all regions the p-values of which exceed the
  threshold \code{cutoff}. By default, this filtering is applied
  to raw p-values (metadata column \dQuote{p.value}). The \code{filterBy}
  argument allows for performing filtering on any of the other three
  p-value metadata columns (if available, otherwise the method quits
  with an error).

  If called for a \code{\linkS4class{GRanges}} object as first argument,
  this method checks if the first argument \code{object} has a metadata
  column named \dQuote{weight.contribution}. If it exists, it returns a
  \code{\linkS4class{GRanges}} object with the elements of \code{object}
  that have a value greater than \code{cutoff} in the
  \dQuote{weight.contribution} metadata column. If this metadata column
  does not exist, the method quits with an error.

  If called for a \code{\linkS4class{GRangesList}} object as first
  argument \code{object}, this method applies the \code{filterResult}
  method for each of its list components and returns a
  \code{\linkS4class{GRangesList}} object. If any of the components
  of \code{object} does not have a metadata column named
  \dQuote{weight.contribution}, the method quits with an error.
}
\value{
  an object of class \code{\linkS4class{AssocTestResultRanges}},
  \code{\linkS4class{GRanges}}, or \code{\linkS4class{GRangesList}} (see
  details above)
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\linkS4class{AssocTestResultRanges}},
  \code{\link[podkat:p.adjust]{p.adjust}}}
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
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
nm.lin <- nullModel(y ~ ., pheno)

## perform association test for multiple regions
res <- assocTest(Z, nm.lin, windows)
res.adj <- p.adjust(res, method="BH")

## show filtered results
res.f <- filterResult(res.adj)
print(res.f)
res.f <- filterResult(res.adj, filterBy="p.value.adj")
print(res.f)

## compute contributions
contrib <- weights(res.f, Z, nm.lin)
contrib

## extract most indicative variants
filterResult(contrib[[1]])
filterResult(contrib)
}
\keyword{methods}
