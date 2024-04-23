\name{sort-methods}
\docType{methods}
\alias{sort}
\alias{sort-methods}
\alias{sort,AssocTestResultRanges-method}
\title{Sort Association Test Results}
\description{
  Rearrange association test results according to sort criterion
}
\usage{
\S4method{sort}{AssocTestResultRanges}(x, decreasing=FALSE,
     sortBy=c("p.value", "p.value.adj", "p.value.resampled",
              "p.value.resampled.adj", "genome"))
}
\arguments{
  \item{x}{object of class \code{\linkS4class{AssocTestResultRanges}}}
  \item{decreasing}{logical indicating if sorting should be done in
    decreasing order}
  \item{sortBy}{sort criterion (see details below)}
}
\details{
  The function \code{sort} takes an \code{\linkS4class{AssocTestResultRanges}}
  object \code{x} and returns a new object of the same class, but with
  the regions rearranged according to the sort criterion \code{sortBy}.
  As an example, if \code{sortBy} is \dQuote{p.value} regions are sorted
  according to raw p-values in ascending order, analogously for the choices
  \dQuote{p.value.adj}, \dQuote{p.value.resampled}, and
  \dQuote{p.value.resampled.adj}. If \code{sortBy} is \dQuote{genome}, the
  regions are arranged along the genome in the same way as a
  \code{\linkS4class{GRanges}} object would be sorted.
  If \code{decreasing} is \code{TRUE}, the order is reversed.
}
\value{
  an \code{\linkS4class{AssocTestResultRanges}} object;
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\linkS4class{AssocTestResultRanges}}}
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

## perform multiple testing correction
res.adj <- p.adjust(res)

## show sorted results (default: by raw p-value)
as(sort(res), "GRanges")
print(sort(res), sortBy="none")

## show results sorted by adjusted p-value
as(sort(res.adj, sortBy="p.value.adj"), "GRanges")
print(sort(res.adj, sortBy="p.value.adj"), sortBy="none")
}
\keyword{methods}
