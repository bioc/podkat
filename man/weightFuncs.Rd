\name{weightFuncs}
\alias{weightFuncs}
\alias{betaWeights}
\alias{logisticWeights}
\alias{invSdWeights}
\title{Weighting Functions}
\description{
  Functions for computing SNV weights from minor allele frequences (MAF)
}
\usage{
betaWeights(x, shape1=1, shape2=25)
logisticWeights(x, th=0.07, slope=150)
invSdWeights(x)
}
\arguments{
  \item{x}{a numeric vector of minor allele frequencies (MAFs); see
    details below}
  \item{shape1,shape2}{shape parameters of Beta distribution weighting
    function (see \code{\link{dbeta}} for details)}
  \item{th,slope}{parameters of the logistic weighting function
    (see details below)}
}
\details{The function \code{betaWeights} is a wrapper around the
  \code{\link{dbeta}} function. It uses the same parameters
  \code{shape1} and \code{shape2}, but does not support the
  non-centrality parameter \code{ncp}. The
  defaults are \code{shape1=1} and \code{shape2=25} as suggested by
  \cite{Wu et al. (2011)} and implemented in the \pkg{SKAT} package.
  If called without argument \code{x}, a function with a
  single argument \code{x} is returned that can directly be used
  as weighting function, e.g. passed as \code{weightFunc} argument to
  the \code{\link{assocTest}} method.

  The function \code{logisticWeights} provides a logistic weighting that
  corresponds to a soft threshold function.
  The \code{th} parameter corresponds to the threshold and
  the \code{slope} parameter corresponds to the steepness of the
  soft threshold.
  Like \code{betaWeights}, this function can be called without \code{x}
  argument to produce a parameter-free weighting function.

  The function \code{invSdWeights} computes weights as
  suggested by \cite{Madsen and Browning (2009)}. For consistency, this function
  also returns a single-argument function if called without \code{x}
  argument.

  For mathematical details, see Subsection 9.3 of the package vignette.
}    
\value{
  a numeric vector with weights as long as the argument \code{x},
  a function if \code{x} was missing;
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}

Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011)
Rare-variant association testing for sequencing data with the sequence
kernel association test. \emph{Am. J. Hum. Genet.} \bold{89}, 82-93.
DOI: \doi{10.1016/j.ajhg.2011.05.029}.

Madsen, B. E. and Browning, S. R. (2009)
A groupwise association test for rare mutations using a weighted
sum statistic. \emph{PLoS Genetics} \bold{5}, e1000384.
DOI: \doi{10.1371/journal.pgen.1000384}}
\seealso{\code{\linkS4class{GenotypeMatrix}}, \code{\link{dbeta}},
\code{\link{assocTest}}}
\examples{
## create a toy example
A <- matrix(rbinom(50, 2, prob=0.2), 5, 10)
MAF <- colSums(A) / (2 * nrow(A))

## compute some weight vectors
betaWeights(MAF, 1, 25)
betaWeights(MAF, 1, 30)
logisticWeights(MAF)
invSdWeights(MAF)

## plot weighting functions (note the missing 'x' arguments)
plot(betaWeights(shape2=30), xlim=c(0, 1))
plot(logisticWeights(), xlim=c(0, 1))
plot(invSdWeights, xlim=c(0, 1))
}
