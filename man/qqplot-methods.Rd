\name{qqplot}
\docType{methods}
\alias{qqplot}
\alias{qqplot-methods}
\alias{qqplot,AssocTestResultRanges,missing-method}
\alias{qqplot,AssocTestResultRanges,AssocTestResultRanges-method}
\title{Quantile-Quantile Plots}
\description{
  Functions for visualizing association test results by means of a
  quantile-quantile (Q-Q) plot
}
\usage{
\S4method{qqplot}{AssocTestResultRanges,missing}(x, y,
       xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
       common.scale=TRUE, preserveLabels=FALSE, lwd=1,
       lcol="red", ..., conf.level=NULL, conf.args=NULL)
\S4method{qqplot}{AssocTestResultRanges,AssocTestResultRanges}(x, y,
       xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
       common.scale=TRUE, preserveLabels=FALSE, lwd=1,
       lcol="red", ..., conf.level=NULL, conf.args=NULL)
}
\arguments{
  \item{x,y}{objects of class
    \code{\linkS4class{AssocTestResultRanges}}}
  \item{xlab}{if \code{preserveLabels} is \code{TRUE}, \code{xlab} is
    interpreted as axis label for the horizontal axis; if
    \code{preserveLabels} is \code{FALSE}, \code{xlab} can be a
    character string or expression that is interpreted as a name/label
    for the object \code{x} and is used for determining an appropriate
    axis label.}
  \item{ylab}{if \code{preserveLabels} is \code{TRUE}, \code{ylab} is
    interpreted as axis label for the vertical axis; if
    \code{preserveLabels} is \code{FALSE}, \code{ylab} can be a
    character string or expression that is interpreted as a name/label
    for the object \code{y} and is used for determining an appropriate
    axis label.}
  \item{common.scale}{if \code{TRUE} (default), the same plotting ranges
    are used for both axes; if \code{FALSE}, the two axes are scaled
    independently.}
  \item{preserveLabels}{if \code{TRUE}, \code{xlab} and \code{ylab} are
    used as axis labels without any change; if \code{FALSE} (default), the
    function interprets \code{xlab} and \code{ylab} as object labels
    for \code{x} and \code{y} and uses them for determining axis labels
    appropriately}
  \item{lwd}{line width for drawing the diagonal line which theoretically
    corresponds to the equality of the two distributions; if zero, no
    diagonal line is drawn.}
  \item{lcol}{color for drawing the diagonal line}
  \item{conf.level}{dummy argument for compatibility with the new
    version of the function in the \pkg{stats} package; currently unused}
  \item{conf.args}{dummy argument for compatibility with the new
    version of the function in the \pkg{stats} package; currently unused}
  \item{...}{all other arguments are passed to
    \code{\link[graphics:plot]{plot}};}
}
\details{
  If \code{qqplot} is called for an
  \code{\linkS4class{AssocTestResultRanges}}
  object without specifying the second argument \code{y},
  a Q-Q plot of the raw p-values in \code{x} against a uniform
  distribution of expected p-values is created, where the theoretical
  p-values are computed using the \code{\link{ppoints}} function.
  In this case, the log-transformed observed p-values contained in \code{x}
  are plotted on the vertical axis and the log-transformed expected
  p-values are plotted
  on the  horizontal axis. If \code{preserveLabels} is \code{TRUE},
  \code{xlab} and \code{ylab} are used as axis labels as usual.
  However, if \code{preserveLabels} is \code{FALSE}, which is the
  default, \code{xlab} is interpreted as object label for \code{x}, i.e.
  the object whose p-values are plotted on the vertical axis.

  If \code{qqplot} is called for two
  \code{\linkS4class{AssocTestResultRanges}} object \code{x} and
  \code{y}, the log-transformed raw p-values of \code{x} and \code{y}
  are plotted against each other, where the p-values of \code{x} are plotted on
  the horizontal axis and the p-values of \code{x} are plotted on the
  vertical axis.
}
\value{
  like the standard \code{\link[stats:qqplot]{qqplot}} function from
  the \pkg{stats} package, \code{qqplot} returns an invisible list
  containing the two sorted vectors of p-values.
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

## perform association tests
res.p <- assocTest(Z, nm.lin, windows, kernel="linear.podkat")
res.s <- assocTest(Z, nm.lin, windows, kernel="linear.SKAT")

## plot results
qqplot(res.p)
qqplot(res.p, res.s, xlab="PODKAT results", ylab="SKAT results")
qqplot(res.p, res.s, xlab="PODKAT results", ylab="SKAT results",
       preserveLabels=TRUE)
qqplot(res.p, res.s, common.scale=FALSE)
}
\keyword{methods}
