\name{AssocTestResultRanges-class}
\docType{class}
\alias{AssocTestResultRanges-class}
\alias{class:AssocTestResultRanges}
\alias{AssocTestResultRanges}
\alias{show,AssocTestResultRanges-method}
\alias{c,AssocTestResultRanges-method}

\title{Class \code{AssocTestResultRanges}}
\description{S4 class for storing the result of an association
  test performed on multiple genomic regions}
\section{Objects}{
  Objects of this class are created by calling \code{\link{assocTest}}
  with a non-empty \code{ranges} argument.
}
\section{Slots}{
This class extends the class \code{\linkS4class{GRanges}} directly and
therefore inherits all its slots and methods.
The following slots are defined for \code{AssocTestResultRanges} objects
additionally:
\describe{
  \item{\code{type}:}{type of null model on which the association test
    was based}
  \item{\code{samples}:}{character vector with sample names (if
    available, otherwise empty)}
  \item{\code{kernel}:}{kernel that was used for the association test}
  \item{\code{weights}:}{weight vector or weighting function that was
    used; \code{NULL} if no weighting was performed}
  \item{\code{width}:}{tolerance radius parameter that was used for
    position-dependent kernels}
  \item{\code{adj.method}:}{which method for multiple testing correction
    has been applied (if any)}
  \item{\code{vcfParams}:}{list of parameters that were used for reading
    genotypes from VCF file}
   \item{\code{sex}:}{factor with sex information (if any)}
 \item{\code{call}:}{the matched call with which the object was created}
}
Apart from these additional slots, all \code{AssocTestResultRanges}
objects have particular metadata columns (accessible via
\code{\link[S4Vectors]{mcols}} or \code{\link[S4Vectors]{elementMetadata}}):
\describe{
  \item{\code{n}:}{number of variants tested in each region; a zero
    does not necessarily mean that there were no variants in this region, it
    only means that no variants were used for testing. Variants are
    omitted from the test if they do not show any variation or if they
    do not satisfy other filter criteria applied by
    \code{\link{assocTest}}. This metadata column is always present.}
  \item{\code{Q}:}{test statistic for each region that was tested. This
    metadata column is always present.}
  \item{\code{p.value}:}{p-value of test for each region that was tested. This
    metadata column is always present.}
  \item{\code{p.value.adj}:}{adjusted p-value of test for each region
    that was tested. This metadata column is only present if multiple
    testing correction has been applied (see \code{\link{p.adjust}}).}
  \item{\code{p.value.resampled}:}{estimated p-value computed as
    the relative frequency of p-values of sampled residuals that
    are at least as significant as the test's p-value in each region.
    This metadata column is only present if resampling has been applied,
    i.e. if \code{\link{assocTest}} has been called with
    \code{n.resampling} greater than zero.}
  \item{\code{p.value.resampled.adj}:}{adjusted empirical p-value (see above).
    This metadata column is only present if resampling and multiple
    testing correction has been applied.}
  }
}
\section{Methods}{
  \describe{
    \item{c}{\code{signature(object="AssocTestResultRanges")}:
      allows for concatenating two or more \code{AssocTestResultRanges}
      objects; this is only meaningful if the different tests have been
      performed on the same samples, on the same genome, with the same
      kernel, and with the same VCF reading parameters (in case that the
      association test has been performed directly on a VCF file).
      All these conditions are checked and if any of them is not
      fulfilled, the method quits with an error. Merging association
      test results that were computed with different \code{sex}
      parameters is possible, but the \code{sex} component is omitted
      and a warning is issued. Note that multiple
      testing correction (see \code{\link{p.adjust}}) should not be
      carried out on parts, but only on the entire set of all tests.
      That is why \code{c} strips off all adjusted p-values.}
    \item{p.adjust}{\code{signature(object="AssocTestResultRanges")}:
      multiple testing correction, see
      \code{\link[podkat:p.adjust-methods]{p.adjust}}.}
    \item{filterResult}{\code{signature(object="AssocTestResultRanges")}:
      apply filtering to p-values or adjusted p-values. For more
      details, see \code{\link{filterResult}}.}
    \item{sort}{\code{signature(object="AssocTestResultRanges")}:
      sort \code{AssocTestResultRanges} object according to specified
      sorting criterion. See \code{\link[podkat:sort-methods]{sort}}
      for more details.}
    \item{plot}{\code{signature(object="AssocTestResultRanges")}:
      make a Manhattan plot of the association test result.
      See \code{\link[podkat:plot-methods]{plot}} for more details.}
    \item{qqplot}{\code{signature(object="AssocTestResultRanges")}:
      make quantile-quantile (Q-Q) plot of association test result.
      See \code{\link[podkat:qqplot-methods]{qqplot}} for more details.}
    \item{show}{\code{signature(object="AssocTestResultRanges")}:
      displays some general information about the result of the
      association test, such as, the number of samples, the number of
      regions tested, the number of regions without variants, the average
      number of variants in the tested regions, the genome, the kernel that
      was applied, and the type of multiple testing correction (if any).}
    \item{print}{\code{signature(x="AssocTestResultRanges")}:
      allows for displaying more information about the object than
      \code{show}. See \code{\link[podkat]{print}}
      for more details.}
  }
}
\section{Accessors and subsetting}{
As mentioned above, the \code{AssocTestResultRanges} inherits all
methods from the \code{\linkS4class{GRanges}} class.
}

\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\link{assocTest}}}
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
print(res.adj)

## show sorted results
as(sort(res.adj), "GRanges")

## show filtered result
print(filterResult(res.adj, cutoff=0.05, filterBy="p.value.adj"))

## make a Manhattan plot
plot(res.adj, which="p.value.adj")
}
\keyword{classes}
