\name{p.adjust-methods}
\docType{methods}
\alias{p.adjust}
\alias{p.adjust-methods}
\alias{p.adjust,AssocTestResultRanges-method}
\title{Adjust p-Value for Multiple Tests}
\description{
  Given an \code{\linkS4class{AssocTestResultRanges}} object,
  this method adds a metadata column with adjusted p-values.
}
\usage{
\S4method{p.adjust}{AssocTestResultRanges}(p, method=p.adjust.methods, n=length(p))
}
\arguments{
  \item{p}{object of class \code{\linkS4class{AssocTestResultRanges}}}
  \item{method}{correction method (see \code{\link{p.adjust.methods}})}
  \item{n}{parameter available for consistency with standard
    \code{\link[stats]{p.adjust}} function; ignored in this implementation}
}
\details{
  This function is a wrapper around the standard
  \code{\link[stats:p.adjust]{p.adjust}} function from the
  \pkg{stats} package. It takes the
  \code{p.value} metadata column from the
  \code{\linkS4class{AssocTestResultRanges}} object \code{p}, applies
  the multiple testing correction method specified as \code{method}
  argument. The method returns a copy of \code{p} with an additional
  metadata column \code{p.value.adj} that contains the adjusted
  p-values. If \code{p} already contained a metadata column
  \code{p.value.adj}, this column is overwritten with the new adjusted
  p-values.

  If \code{p} also contains a metadata column \code{p.value.resampled},
  multiple testing correction is also applied to resampled p-values.
  The resulting adjusted p-values are placed in the metadata column
  \code{p.value.resampled.adj}.

  Note that, for consistency with the standard
  \code{\link[stats:p.adjust]{p.adjust}} function, the default
  correction method is \dQuote{holm}.
}
\value{
  an \code{\linkS4class{AssocTestResultRanges}} object (see details above)
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\linkS4class{AssocTestResultRanges}},
  \code{\link[stats:p.adjust]{p.adjust}}}
\examples{
## load genome description
data(hgA)

## partition genome into overlapping windows
windows <- partitionRegions(hgA)

## load genotype data from VCF file
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
Z <- readGenotypeMatrix(vcfFile)

## read phenotype data from CSV file (continous trait + covariates)
phenoFile <- system.file("examples/example1lin.csv", package="podkat")
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
nm.lin <- nullModel(y ~ ., pheno)

## perform association test for multiple regions
res <- assocTest(Z, nm.lin, windows)
print(res)

## perform multiple testing correction
res.adj <- p.adjust(res, method="BH")
print(res.adj)
}
\keyword{methods}
