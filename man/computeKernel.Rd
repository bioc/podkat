\name{computeKernel}
\alias{computeKernel}
\title{Compute Kernel Matrix}
\description{
  Computes kernel matrix for a given genotype matrix
}
\usage{
computeKernel(Z, kernel=c("linear.podkat", "quadratic.podkat",
              "localsim.podkat", "linear.SKAT", "quadratic.SKAT",
              "localsim.SKAT"), weights=NULL, pos=NULL, width=1000)
}
\arguments{
  \item{Z}{a matrix or an object of class \code{\linkS4class{Matrix}}
    (note that the latter also includes objects of class
    \code{\linkS4class{GenotypeMatrix}})}
  \item{kernel}{type of kernel to use}
  \item{weights}{numeric vector with variant weights;
    must be as long as the number of columns of \code{Z}.
    Use \code{NULL} for unweighted kernels.}
  \item{pos}{numeric vector with positions of variants;
    must be as long as the number of columns of \code{Z}.
    This argument is mandatory for the position-dependent
    kernels \dQuote{linear.podkat}, \dQuote{quadratic.podkat},
    and \dQuote{localsim.podkat}; ignored for kernels
    \dQuote{linear.SKAT}, \dQuote{quadratic.SKAT},
    and \dQuote{localsim.SKAT}.}
  \item{width}{tolerance radius parameter for position-dependent kernels
    \dQuote{linear.podkat}, \dQuote{quadratic.podkat},
    and \dQuote{localsim.podkat} (see details below);
    must be single positive numeric value.
    Ignored for kernels \dQuote{linear.SKAT}, \dQuote{quadratic.SKAT},
    and \dQuote{localsim.SKAT}.}  
}
\details{This function computes a kernel matrix for a given genotype
  matrix \code{Z} and a given kernel. It supposes that \code{Z} is a
  matrix-like object (a numeric matrix, a sparse matrix, or an object of
  class \code{\linkS4class{GenotypeMatrix}}) in which rows correspond
  to samples and columns correspond to variants. There are six different
  kernels available: \dQuote{linear.podkat}, \dQuote{quadratic.podkat},
  \dQuote{localsim.podkat}, \dQuote{linear.SKAT}, \dQuote{quadratic.SKAT},
  and \dQuote{localsim.SKAT}. All of these kernels can be used with or
  without weights. The weights can be specified with the \code{weights}
  argument which must be a numeric vector with as many elements as the
  matrix \code{Z} has columns. If no weighting should be used,
  \code{weights} must be set to \code{NULL}.

  The position-dependent kernels \dQuote{linear.podkat},
  \dQuote{quadratic.podkat}, and \dQuote{localsim.podkat} require the
  positions of the variants in \code{Z}. So, if any of these three
  kernels is selected, the argument \code{pos} is mandatory and must
  be a numeric vector with as many elements as the
  matrix \code{Z} has columns.

  If the \code{pos} argument is \code{NULL} and \code{Z} is a
  \code{\linkS4class{GenotypeMatrix}} object, the positions in
  \code{variantInfo(Z)} are taken. In this case, all variants need to reside
  on the same chromosome. If the variants in \code{variantInfo(Z)} are from
  multiple chromosomes, \code{computeKernel} quits with an error.
  As said, this only happens if \code{pos} is \code{NULL}, otherwise
  the \code{pos} argument has priority over the information stored in
  \code{variantInfo(Z)}.

  For details on how the kernels compute the pairwise similarities of
  genotypes, see Subsection 9.2 of the package vignette.
}
\value{
  a positive semi-definite kernel matrix with as many rows and
  columns as \code{Z} has rows
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}

Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011)
Rare-variant association testing for sequencing data with the sequence
kernel association test. \emph{Am. J. Hum. Genet.} \bold{89}, 82-93.
DOI: \href{http://dx.doi.org/10.1016/j.ajhg.2011.05.029}{10.1016/j.ajhg.2011.05.029}.
}
\seealso{\code{\linkS4class{GenotypeMatrix}}}
\examples{
## create a toy example
A <- matrix(rbinom(50, 2, prob=0.2), 5, 10)
pos <- sort(sample(1:10000, ncol(A)))

## compute some unweighted kernels
computeKernel(A, kernel="linear.podkat", pos=pos, width=100)
computeKernel(A, kernel="localsim.podkat", pos=pos, width=100)
computeKernel(A, kernel="linear.SKAT")

## compute some weighted kernels
MAF <- colSums(A) / (2 * nrow(A))
weights <- betaWeights(MAF)
computeKernel(A, kernel="linear.podkat", pos=pos, weights=weights)
computeKernel(A, kernel="linear.SKAT", weights=weights)
computeKernel(A, kernel="localsim.SKAT", weights=weights)
}
