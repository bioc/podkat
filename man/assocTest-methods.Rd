\name{assocTest}
\docType{methods}
\alias{assocTest}
\alias{assocTest-methods}
\alias{assocTest,GenotypeMatrix,NullModel-method}
\alias{assocTest,matrix,NullModel-method}
\alias{assocTest,TabixFile,NullModel-method}
\alias{assocTest,character,NullModel-method}
\title{Perform Association Test}
\description{
  Method for performing a kernel-based association test given a
  genotype, VCF file, or kernel matrix
}
\usage{
\S4method{assocTest}{GenotypeMatrix,NullModel}(Z, model, ranges,
          kernel=c("linear.podkat", "localsim.podkat",
                   "quadratic.podkat", "linear.SKAT",
                   "localsim.SKAT", "quadratic.SKAT"),
          width=1000, weights=NULL, weightFunc=betaWeights(),
          method=NULL, adj=c("automatic", "none", "force"),
          pValueLimit=0.05)
\S4method{assocTest}{matrix,NullModel}(Z, model, method=NULL,
          adj=c("automatic", "none", "force"), pValueLimit=0.05)
\S4method{assocTest}{TabixFile,NullModel}(Z, model, ranges,
          kernel=c("linear.podkat", "localsim.podkat",
                   "quadratic.podkat", "linear.SKAT",
                   "localsim.SKAT", "quadratic.SKAT"),
          cl=NULL, nnodes=1, batchSize=20,
          noIndels=TRUE, onlyPass=TRUE, na.limit=1, MAF.limit=1,
          na.action=c("impute.major", "omit"),
          MAF.action=c("invert", "omit", "ignore"),
          sex=NULL, weightFunc=betaWeights(), width=1000,
          method=NULL, adj=c("automatic", "none", "force"),
          pValueLimit=(0.1 / length(ranges)), tmpdir=tempdir(),
          displayProgress=TRUE)
\S4method{assocTest}{character,NullModel}(Z, model, ...)
}
\arguments{
  \item{Z}{an object of class \code{\linkS4class{GenotypeMatrix}}, a
    quadratic kernel matrix, an object of class
    \code{\link[Rsamtools]{TabixFile}}, or a character string
    with a file name}
  \item{model}{an object of class \code{\linkS4class{NullModel}}}
  \item{ranges}{an object with genomic regions to be tested; may be
    an object of class \code{\linkS4class{GRanges}} or
    \code{\linkS4class{GRangesList}}. If missing, \code{assocTest} takes the
    whole genotype matrix or the genotypes in the VCF file as a whole.}
  \item{kernel}{determines the kernel that should be used for
    association testing (see Subsection 9.2 of the package vignette for
    details)}
  \item{width}{tolerance radius parameter for position-dependent kernels
    \dQuote{linear.podkat}, \dQuote{quadratic.podkat},
    and \dQuote{localsim.podkat}; must be single positive numeric value;
    ignored for kernels \dQuote{linear.SKAT}, \dQuote{quadratic.SKAT},
    and \dQuote{localsim.SKAT}  (see Subsection 9.2 of the package
    vignette for details)}
  \item{weights}{for the method with signature
    \code{GenotypeMatrix,NullModel}, it is also possible to supply
    weights directly as a numeric vector that is as long as the number of
    columns of \code{Z}. In this case, the argument \code{weightFunc}
    is ignored. Use \code{NULL} (default) to use automatic weighting
    with the function supplied as argument \code{weightFunc}. If
    \code{weightFunc} is \code{NULL} too, no weighting takes place, i.e.
    an unweighted kernel is used.}
  \item{weightFunc}{function for computing variant weights from minor allele
    frequencies (MAFs); see \code{\link{weightFuncs}} for weighting
    and Subsection 9.3 of the package vignette for functions provided by
    the \pkg{podkat} package. Use \code{NULL}
    for unweighted kernels.}
  \item{method}{identifies the method for computing the p-values. If the
    null model is of type \dQuote{logistic} and small sample correction
    is applied (see argument \code{adj} below), possible values are
    \dQuote{unbiased}, \dQuote{population}, \dQuote{sample}, and
    \dQuote{SKAT} (see details below and Subsection 9.5 of the package
    vignette). If the null model is of type \dQuote{linear} or if
    the null model is of type \dQuote{logistic} and no small sample
    correction is applied, possible values are \dQuote{davies},
    \dQuote{liu}, and \dQuote{liu.mod} (see details
    below and Subsection 9.1 of the package vignette). If the null
    model is of type \dQuote{bernoulli}, this argument is ignored.}
  \item{adj}{whether or not to use small sample correction for
    logistic models (binary trait with covariates). The choice
    \dQuote{none} turns off small sample correction. If \dQuote{force}
    is chosen, small sample correction is turned on unconditionally.
    If \dQuote{automatic} is chosen (default), small sample correction
    is turned on if the number of samples does not exceed 2,000. This
    argument is ignored for any type of model except \dQuote{logistic}
    and small sample correction is switched off. For details how to
    train a null model for small sample correction, see
    \code{\link{nullModel}} and Sections 4 and 9.5 of the package
    vignette. An adjustment of higher moments is
    performed whenever sampled null model residuals are available
    in the null model \code{model} (slot \code{res.resampled.adj}, see
    \code{\linkS4class{NullModel}}).}
  \item{pValueLimit}{if the null model is of type \dQuote{bernoulli},
    \code{assocTest} performs an exact mixture of Bernoulli test. This
    test uses a combinatorial algorithm to compute exact p-values and,
    for the sake of computational efficiency, quits if a pre-specified
    p-value threshold is exceeded. This threshold can be specified with
    the \code{pValueLimit} argument. This argument is ignored for other
    types of tests/null models.}
  \item{cl}{if \code{cl} is an object of class
    \code{\link[parallel:makePSOCKcluster]{SOCKcluster}}, association
    testing is carried out in parallel on the cluster specified by \code{cl}.
    If \code{NULL} (default), either no parallelization is done (if
    \code{nnodes=1}) or \code{assocTest} launches a cluster with
    \code{nnodes} R client processes on \code{localhost}. See Subsection
    8.5.2 of the package vignette.}
  \item{nnodes}{if \code{cl} is \code{NULL} and \code{nnodes} is greater
    than 1, \code{\link{makePSOCKcluster}} is called with \code{nnodes}
    nodes on \code{localhost}, i.e. \code{nnodes} R slave processes are
    launched on which association testing is carried out in parallel.
    The default is \code{1}. See Subsection 8.5.2 of
    the package vignette.}
  \item{batchSize}{parameter which determines how many
    regions of \code{ranges}
    are processed at once. The larger \code{batchSize}, the larger the
    the batches that are read from the VCF file \code{Z}. A larger
    \code{batchSize} reduces the number of individual read operations,
    which improves performance. However, a larger \code{batchSize} also
    requires larger amounts of memory. A good choice of \code{batchSize},
    therefore, depends on the size and sparseness of the VCF file and as
    well on the available memory. See Subsection 8.5 of
    the package vignette.}
  \item{noIndels}{if \code{TRUE} (default), only single nucleotide
    variants (SNVs) are considered
    and indels in the VCF file \code{Z} are skipped.}
  \item{onlyPass}{if \code{TRUE} (default), only variants are considered
    whose value in the \code{FILTER} column is \dQuote{PASS}.}
   \item{na.limit}{all variants with a missing value ratio above this
    threshold in the VCF file \code{Z} are not considered.}
  \item{MAF.limit}{all variants with a minor allele frequency (MAF) above this
    threshold in the VCF file \code{Z} are not considered.}
  \item{na.action}{if \dQuote{impute.major}, all missing values will
    be imputed by major alleles before association testing. If
    \dQuote{omit}, all columns containing missing values in the VCF file
    \code{Z} are ignored.}
  \item{MAF.action}{if \dQuote{invert}, all columns with an MAF
    exceeding 0.5 will be inverted in the sense that all minor alleles
    will be replaced by major alleles and vice versa.
    If \dQuote{omit}, all variants in the VCF file with an MAF greater
    than 0.5 are ignored. If \dQuote{ignore},
    no action is taken and MAFs greater than 0.5 are kept as they are.}
  \item{sex}{if \code{NULL}, all samples are treated the same
    without any modifications; if \code{sex} is a factor with levels
    \code{F} (female) and \code{M} (male) that is as long as
    \code{length(model)}, this argument is
    interpreted as the sex of the samples. In this case, the genotypes
    corresponding to male samples are
    doubled before further processing. This is designed for mixed-sex
    analyses of the X chromosome outside of the pseudoautosomal
    regions.}
  \item{tmpdir}{if computations are parallelized over multiple client
    processes (see arguments \code{nnodes} and \code{cl}), the exchange
    of the null model object between the master process and the client
    processes is done via a temporary file. The \code{tmpdir} argument
    allows to specify into which directory the temporary file should be
    saved. On multi-core systems, the default should be sufficient.
    If the computations are distributed over a custom cluster,
    the \code{tmpdir} argument needs to be chosen such that all clients
    can access it via the same path.}
  \item{displayProgress}{if \code{TRUE} (default) and if \code{ranges} is a
    \code{\linkS4class{GRangesList}}, a progress message is printed
    upon completion of each list component (typically consisting of
    regions of one chromosome); this argument is ignored if
    \code{ranges} is not an object of class \code{\linkS4class{GRangesList}}.}
  \item{...}{all other parameters are passed on to the \code{assocTest}
    method with signature \code{TabixFile,NullModel}.}
}
\details{
  The \code{assocTest} method is the main function of the \pkg{podkat}
  package. For a given genotype and a null model, it performs the actual
  association test(s).

  For null models of types \dQuote{linear} and
  \dQuote{logistic} (see \code{\linkS4class{NullModel}} and
  \code{\link{nullModel}}), a \emph{variance component score test} is used
  (see Subsection 9.1 of the package vignette for details).
  The test relies on the choice of a particular kernel to measure the
  pairwise similarities of genotypes. The choice of the kernel can be
  made with the \code{kernel} argument (see \code{\link{computeKernel}}
  and Subsection 9.2 of the package vignette for more details).
  For null models of  type \dQuote{linear}, the test statistic follows
  a \emph{mixture of chi-squares distribution}. For models of typ
  \dQuote{logistic}, the test statistic approximately follows a mixture
  of chi-squares distribution. The computation of p-values for a given
  mixture of chi-squares can be done according to \cite{Davies (1980)} (which is
  the default), according to \cite{Liu et al. (2009)}, or using a
  modified method similar to the one suggested by
  \cite{Liu et al. (2009)} as implemented in the \pkg{SKAT} package,
  too. Which method is used can be controlled using the \code{method}
  argument. If method according to \cite{Davies (1980)} fails,
  \code{assocTest} resorts to the method by \cite{Liu et al. (2009)}.
  See also Subsection 9.1 of the package vignette for more details.

  For null models of type \dQuote{logistic}, the \code{assocTest}
  method also offers the small sample correction suggested by
  \cite{Lee et al. (2012)}. Whether small sample correction is applied,
  is controlled by the \code{adj} argument. The additional adjustment
  of higher moments as suggested by \cite{Lee et al. (2012)} is
  performed whenever resampled null model residuals are available in
  the null model \code{model} (slot \code{res.resampled.adj}, see
  \code{\linkS4class{NullModel}}). In this case, the \code{method}
  argument controls how the excess kurtosis of test statistics sampled
  from the null distribution are computed. The
  default setting \dQuote{unbiased} computes unbiased estimates by
  using the exact expected value determined from the mixture components.
  The settings \dQuote{population} and
  \dQuote{sample} use almost unbiased and biased sample statistics,
  respectively. The choice \dQuote{SKAT} uses the same method as
  implemented in the \pkg{SKAT} package. See Subsection 9.5 of the
  package vignette for more details.
  
  If the null model is of type \dQuote{bernoulli}, the test statistic
  follows a mixture of Bernoulli distributions. In this case, an exact
  p-value is determined that is computed as the probability to observe
  a test statistic for random Bernoulli-distributed traits (under the null
  hypothesis) that is at least as large as the observed test
  statistic. For reasons of computational complexity, this option is
  limited to sample numbers not larger than 100. See Subsection 9.1 of
  the package vignette for more details.

  The \pkg{podkat} package offers multiple interfaces for association
  testing all of which require the second argument \code{model} to be a
  \code{\linkS4class{NullModel}} object.
  The simplest method is to call \code{assocTest} for an object
  of class \code{\linkS4class{GenotypeMatrix}} as first argument
  \code{Z}. If the \code{ranges} argument is not supplied, a single
  association test is performed using the entire genotype matrix
  contained in \code{Z} and an object of class
  \code{\linkS4class{AssocTestResult}} is returned. In this case,
  all variants need to reside on the same chromosome (compare with
  \code{\link{computeKernel}}). If the \code{ranges} argument is
  specified, each region in \code{ranges} is tested separately and the
  result is returned as an \code{\linkS4class{AssocTestResultRanges}}
  object.

  As said, the simplest method is to store the entire genotype in a
  \code{\linkS4class{GenotypeMatrix}} object and to call \code{assocTest}
  as described above. This approach has the shortcoming that the entire
  genotype must be read (e.g. from a VCF file) and kept in memory as a
  whole. For large studies, in particular, whole-genome studies, this is
  not feasible. In order to be able to cope with large studies, the
  \pkg{podkat} package offers an interface that allows for reading from
  a VCF file piece by piece without the need to read and store the entire
  genotype at once. If \code{Z} is a \code{\link[Rsamtools]{TabixFile}}
  object or the name of a VCF file, \code{assocTest} reads from the file
  in batches of \code{batchSize} regions, performs the association tests
  for these regions, and returns the results as an
  \code{\linkS4class{AssocTestResultRanges}} object. This sequential
  batch processing can also be parallelized. The user can either set up
  a cluster him-/herself and pass the 
  \code{\link[parallel:makePSOCKcluster]{SOCKcluster}} object as
  \code{cl} argument. If the \code{cl} is \code{NULL}, users can
  leave the setup of the cluster to \code{assocTest}. In this case, the
  only thing necessary is to determine the number of R client processes
  by the \code{nnodes} argument. The variant with the VCF interface
  supports the same pre-processing and filter arguments as
  \code{\link{readGenotypeMatrix}} to control which variants are actually
  taken into account and how to handle variants with MAFs greater than
  50\%.
  
  If the argument \code{Z} is a numeric matrix, \code{Z} is
  interpreted as a kernel matrix
  \eqn{\boldmath{K}}{\bold{K}}. Then a single association test is
  performed as described above and the result is returned as an
  \code{\linkS4class{AssocTestResult}} object. This allows the user to
  use a custom kernel not currently implemented in the \pkg{podkat}
  package. The \code{assocTest} function assumes that row and column
  objects in the kernel matrix are in the same order. It does not
  perform any check whether row and column names are the same or whether
  the kernel matrix is actually positive semi-definite. Users should be
  aware that running the function for invalid kernels matrices,
  i.e. for a matrix that is not positive semi-definite, produces
  meaningless results and may even lead to unexpected errors.

  Finally, note that the samples in the null model \code{model} and
  in the genotype (\code{GenotypeMatrix} object or VCF file)
  need not be aligned to each other. If both the samples in \code{model}
  and in the genotype are named (i.e. row names are defined
  for \code{Z} if it is a \code{GenotypeMatrix} object; VCF files
  always contain sample names anyway), \code{assocTest} checks if all
  samples in \code{model} are present in the genotype. If so, it
  selects only those samples from the genotype that occur in the null
  model. If not, it quits with an error. If either the samples in the
  null model or the genotypes are not named, \code{assocTest} assumes
  that the samples are aligned to each other. This applies only if
  the number of samples in the null model and the number of genotypes
  are the same or if the number of genotypes equals the number of
  samples in the null model plus the number of samples that were omitted
  from the null model when it was trained (see
  \code{\linkS4class{NullModel}} and \code{\link{nullModel}}).
  Otherwise, the function quits with an error.
  An analogous procedure is applied if the kernel matrix interface is
  used (signature \code{matrix,NullModel}).
}
\value{
  an object of class \code{\linkS4class{AssocTestResult}} or
  \code{\linkS4class{AssocTestResultRanges}} (see details above)
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}

Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011)
Rare-variant association testing for sequencing data with the sequence
kernel association test. \emph{Am. J. Hum. Genet.} \bold{89}, 82-93.
DOI:
\href{http://dx.doi.org/10.1016/j.ajhg.2011.05.029}{10.1016/j.ajhg.2011.05.029}.

Lee, S., Emond, M. J., Bamshad, M. J., Barnes, K. C., Rieder, M. J.,
Nickerson, D. A., NHLBI Exome Sequencing Project - ESP Lung Project
Team, Christiani, D. C., Wurfel, M. M., and Lin, X. (2012)
Optimal unified approach for rare-variant association testing with
application to small-sample case-control whole-exome sequencing studies.
\emph{Am. J. Hum. Genet.} \bold{91}, 224-237. DOI:
\href{http://dx.doi.org/10.1016/j.ajhg.2012.06.007}{10.1016/j.ajhg.2012.06.007}.

Davies, R. B. (1980) The distribution of a linear combination of
\eqn{\chi^2}{chi^2} random
variables. \emph{J. R. Stat. Soc. Ser. C-Appl. Stat.} \bold{29},
323-333.

Liu, H., Tang, Y., and Zhang, H. (2009) A new chi-square
approximation to the distribution of non-negative definite quadratic
forms in non-central normal variables. \emph{Comput. Stat. Data Anal.}
\bold{53}, 853-856.
}
\seealso{\code{\link{AssocTestResult}},
  \code{\link{AssocTestResultRanges}}, \code{\link{nullModel}},
  \code{\linkS4class{NullModel}}, \code{\link{computeKernel}},
  \code{\link{weightFuncs}}, \code{\link{readGenotypeMatrix}},
  \code{\linkS4class{GenotypeMatrix}},
  \code{\link{plot}}, \code{\link{qqplot}}, \code{\link{p.adjust}},
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
pheno.c <- read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
model.c <- nullModel(y ~ ., pheno.c)

## perform association test
res <- assocTest(Z, model.c, windows)
print(res)

## perform association test using the VCF interface
res <- assocTest(vcfFile, model.c, windows, batchSize=100)
print(res)

## create Manhattan plot of adjusted p-values
plot(p.adjust(res), which="p.value.adj")

## read phenotype data from CSV file (binary trait + covariates)
phenoFile <- system.file("examples/example1log.csv", package="podkat")
pheno.b <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
model.b <- nullModel(y ~ ., pheno.b)

## perform association test
res <- assocTest(Z, model.b, windows)
print(res)

## create Manhattan plot of adjusted p-values
plot(p.adjust(res), which="p.value.adj")
}
\keyword{methods}
