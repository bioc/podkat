\name{print-methods}
\docType{methods}
\alias{print}
\alias{print-methods}
\alias{print,AssocTestResultRanges-method}
\title{Print Association Test Results}
\description{
  Display method for S4 class \code{\linkS4class{AssocTestResultRanges}}
}
\usage{
\S4method{print}{AssocTestResultRanges}(x, cutoff=0.05,
        sortBy=c("p.value", "p.value.adj", "p.value.resampled",
                 "p.value.resampled.adj", "genome", "none"),
        max.show=10)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{AssocTestResultRanges}}}
  \item{cutoff}{a numerical vector with one or more p-value thresholds;
    if present (otherwise \code{NA} or an empty vector must be passed),
    \code{print} displays the number of tested regions with a p-value
    below each threshold. If the
    \code{\linkS4class{AssocTestResultRanges}} object also contains
    adjusted p-values, the numbers of tested regions with p-values
    below each of the thresholds are printed too. If \code{max.show} is
    greater than 0, the \code{max.show} most significant regions up to
    an (adjusted) p-value (depending on the \code{sortBy} argument)
    up to the largest threshold are shown.}
  \item{sortBy}{a character string that determines (1) how regions are
    sorted and (2) according to which p-value the \code{cutoff}
    threshold is applied when printing regions; if \code{sortBy} is
    \dQuote{p.value} (default), regions are sorted according to
    raw p-values and only the \code{max.show} most significant
    regions are printed - as long as the raw p-value is not larger than
    the largest value in the \code{cutoff} argument. For
    \dQuote{p.value.adj}, regions are sorted and filtered
    according to adjusted p-values, analogously for choices
    \dQuote{p.value.resampled} and
    \dQuote{p.value.resampled.adj}. In case that \code{sortBy}
    is \dQuote{genome}, the p-values are ignored and
    the first \code{max.show} regions in the genome are displayed.
    In case that \code{sortBy}
    is \dQuote{none}, the p-values are also ignored and
    the first \code{max.show} regions are displayed in the order as they
    appear in the \code{\linkS4class{AssocTestResultRanges}} object.
  }
  \item{max.show}{maximum number of regions to display; if 0, no regions
    are displayed at all.}
}
\details{
  \code{print} displays the most important information stored in
  an \code{\linkS4class{AssocTestResultRanges}} object \code{x}.
  That includes the   type of null model, the numbers of samples and
  tested regions, the kernel that was used for testing, etc. Depending
  on the \code{cutoff} argument, a certain number of significant tests
  is printed. If \code{max.show}
  is larger than 0, then some regions are shown along with association
  test results. Which regions are selected and how they are sorted,
  depends on the arguments \code{sortBy} and \code{cutoff} (see above).
}
\value{
  \code{print} returns its argument \code{x} invisibly.
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\linkS4class{GenotypeMatrix}},
  \code{\linkS4class{NullModel}}, \code{\linkS4class{AssocTestResult}},
  \code{\link{AssocTestResultRanges}}
}
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

## print 'AssocTestResultRanges' object
print(res)
print(res, max.show=0)
print(res.adj, cutoff=c(0.005, 0.01, 0.05))
print(res.adj, cutoff=0.05, sortBy="p.value.adj")
print(res.adj, cutoff=NA, sortBy="none", max.show=40)
}
\keyword{methods}
