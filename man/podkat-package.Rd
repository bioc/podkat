\name{podkat-package}
\alias{podkat-package}
\alias{podkat}
\docType{package}
\title{PODKAT Package}
\description{
  This package provides an association test that is capable of
  dealing with very rare and even private variants.
  This is accomplished by a kernel-based approach that takes
  the positions of the variants into account. The test can be used
  for pre-processed matrix data, but also directly for variant
  data stored in VCF files. Association testing can be performed
  whole-genome, whole-exome, or restricted to pre-defined regions
  of interest. The test is complemented by tools for analyzing
  and visualizing the results.}
\details{The central method of this package is \code{\link{assocTest}}.
  It provides several different kernel-based association tests, in
  particular, the position-dependent kernel association test (PODKAT),
  but also some variants of the SNP-set kernel association test (SKAT).
  The test can be run for genotype data given in (sparse) matrix format
  as well as directly on genotype data stored in a variant call format
  (VCF) file. In any case, the user has to create a null model by
  the \code{\link{nullModel}} function beforehand. Upon completion of
  an association test, the package also provides methods for filtering,
  sorting, multiple testing correction, and visualization of results.
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}}
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

## perform association test
res <- assocTest(Z, nm.lin, windows)

## display results
print(res)
print(p.adjust(res))
plot(p.adjust(res), which="p.value.adj")
}
\keyword{package}
