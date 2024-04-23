\name{AssocTestResult-class}
\docType{class}
\alias{AssocTestResult-class}
\alias{class:AssocTestResult}
\alias{AssocTestResult}
\alias{show,AssocTestResult-method}

\title{Class \code{AssocTestResult}}
\description{S4 class for storing the result of an association
  test for a single genomic region}
\section{Objects}{
  Objects of this class are created by calling \code{\link{assocTest}}
  for a single genomic region.
}
\section{Slots}{
The following slots are defined for \code{AssocTestResult} objects:
\describe{
  \item{\code{type}:}{type of null model on which the association test
    was based}
  \item{\code{samples}:}{character vector with sample names (if
    available, otherwise empty)}
  \item{\code{kernel}:}{kernel that was used for the association test}
  \item{\code{dim}:}{dimensions of genotype matrix that was tested}
  \item{\code{weights}:}{weight vector that was used; empty if no
    weighting was performed}
  \item{\code{width}:}{tolerance radius parameter that was used for
    position-dependent kernels}
  \item{\code{method}:}{method(s) used to compute p-values; a single
    character string if no resampling was done, otherwise a list with
    two components specifying the p-value computation method for the
    test's p-value and the resampled p-values separately.}
  \item{\code{correction}:}{a logical vector indicating whether the
    small sample correction was carried out (first component
    \code{exact} is \code{TRUE}) and/or higher moment correction was carried
    out (second component \code{resampling} is \code{TRUE}).}
  \item{\code{Q}:}{test statistic}
  \item{\code{p.value}:}{the test's p-value}
  \item{\code{Q.resampling}:}{test statistics for sampled null model
    residuals}
  \item{\code{p.value.resampling}:}{p-values for sampled null model
    residuals}
  \item{\code{p.value.resampled}:}{estimated p-value computed as
    the relative frequency of p-values of sampled residuals that
    are at least as significant as the test's p-value}
  \item{\code{call}:}{the matched call with which the object was created}
}}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="AssocTestResult")}:
      displays the test statistic and the p-value along with the type of
      the null model, the number of samples, the number of SNVs, and the
      kernel that was used to carry out the test.}
  }
}

\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
\seealso{\code{\link{assocTest}}}
\examples{
## load genome description
data(hgA)

## load genotype data from VCF file
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
Z <- readGenotypeMatrix(vcfFile)

## read phenotype data from CSV file (continuous trait + covariates)
phenoFile <- system.file("examples/example1lin.csv", package="podkat")
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
nm.lin <- nullModel(y ~ ., pheno)

## perform association test for entire genotype matrix
res <- assocTest(Z, nm.lin)
show(res)

## perform association test for subset of genotype matrix
res <- assocTest(Z[, 50:100], nm.lin)
show(res)
}
\keyword{classes}
