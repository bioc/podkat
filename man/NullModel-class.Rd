\name{NullModel-class}
\docType{class}
\alias{NullModel-class}
\alias{class:NullModel}
\alias{NullModel}
\alias{residuals,NullModel-method}
\alias{names,NullModel-method}
\alias{show,NullModel-method}
\alias{length,NullModel-method}
\alias{coefficients,NullModel-method}
\alias{[,NullModel,index,missing,missing-method}

\title{Class \code{NullModel}}
\description{S4 class for storing null models for later usage with the
  \code{\link{assocTest}} method}
\section{Objects}{
  Objects of this class are created by calling \code{\link{nullModel}}.
}
\section{Slots}{
The following slots are defined for \code{NullModel} objects:
    \describe{
    \item{\code{type}:}{type of model}
    \item{\code{residuals}:}{residuals of linear model; for type
      \dQuote{bernoulli}, this is simply the trait vector (see
      \code{\link{nullModel-methods}} for details)}
    \item{\code{model.matrix}:}{model matrix of the (generalized) linear
      model trained for the covariates (if any)}
    \item{\code{inv.matrix}:}{pre-computed inverse of some matrix
      needed for computing the null distribution; only used for
      types \dQuote{logistic} and \dQuote{linear}}
    \item{\code{P0sqrt}:}{pre-computed square root of matrix \eqn{P_0}
      (see Subsections 9.1 and 9.5 of the package vignette);
      needed for computing the null distribution in case the small
      sample correction is used for a logistic model; computed only
      if \code{\link{nullModel}} is called with \code{adjExact=TRUE}.}
    \item{\code{coefficients}:}{coefficients of (generalized) linear
      model trained for the covariates (if any)}
    \item{\code{na.omit}:}{indices of samples omitted from (generalized)
      linear model because of missing values in target or covariates}
    \item{\code{n.cases}:}{for binary traits (types \dQuote{logistic}
      and \dQuote{bernoulli}), the number of cases, i.e. the number of 1's
      in the trait vector}
    \item{\code{variance}:}{for continuous traits (type
      \dQuote{linear}), this is a single numeric value with the variance
      of residuals of the linear model; for logistic models with
      binary traits (type \dQuote{logistic}), this is a vector with
      variances of the per-sample Bernoulli distributions; for later
      use of the exact mixture-of-Bernoulli test (type \dQuote{bernoulli}),
      this is the variance of the Bernoulli distribution}
    \item{\code{prob}:}{for logistic models with
      binary traits (type \dQuote{logistic}), this is a vector with
      probabilities of the per-sample Bernoulli distributions; for later
      use of the exact mixture-of-Bernoulli test (type \dQuote{bernoulli}),
      this is the probability of the Bernoulli distribution}
     \item{\code{type.resampling}:}{which resampling algorithm was used}
     \item{\code{res.resampling}:}{matrix with residuals sampled under
       the null hypothesis (if any)}
    \item{\code{res.resampling.adj}:}{matrix with residuals sampled under
       the null hypothesis for the purpose of higher moment correction (if
       any; only used for logistic models with small sample correction)}
   \item{\code{call}:}{the matched call with which the object was created}
  }
}
\section{Details}{
  This class serves as the general interface for storing the necessary
  phenotype information for a later association test. Objects of this
  class should only be created by the \code{\link{nullModel}} function.
  Direct modification of object slots is strongly discouraged!
}
\section{Methods}{
   \describe{
     \item{show}{\code{signature(object="NullModel")}:
       displays basic information about the null model, such as,
       the type of the model and the numbers of covariates.}
     }
}
\section{Accessors}{
   \describe{
     \item{residuals}{\code{signature(object="NullModel")}:
       returns the \code{residuals} slot.}
     \item{names}{\code{signature(object="NullModel")}:
       returns the names of samples in the null model.}
     \item{coefficients}{\code{signature(object="NullModel")}:
       returns the \code{coefficients} slot.}
     \item{length}{\code{signature(x="NullModel")}:
       returns the number of samples that was used to train the null model.}
    }
}
\section{Subsetting}{
  For a \code{NullModel} object \code{x} and an index vector \code{i}
  that is a permutation of \code{1:length(x)},
  \code{x[i]} returns a new \code{NullModel} object in which the samples
  have been rearranged according to the permutation \code{i}. This is
  meant for applications in which the order of the samples in a subsequent
  association test is different from the order of the samples when the
  null model was trained/created.
}

\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\link{nullModel}}}
\examples{
## read phenotype data from CSV file (continuous trait + covariates)
phenoFile <- system.file("examples/example1lin.csv", package="podkat")
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
model <- nullModel(y ~ ., pheno)
model
length(model)
residuals(model)

## read phenotype data from CSV file (binary trait + covariates)
phenoFile <- system.file("examples/example1log.csv", package="podkat")
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
model <- nullModel(y ~ ., pheno)
model
length(model)
residuals(model)

## "train" simple Bernoulli model on a subset of 100 samples
model <- nullModel(y ~ 0, pheno[1:100, ])
model
length(model)
residuals(model)
}
\keyword{classes}
