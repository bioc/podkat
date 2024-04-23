\name{nullModel}
\docType{methods}
\alias{nullModel}
\alias{nullModel-methods}
\alias{nullModel,matrix,numeric-method}
\alias{nullModel,matrix,factor-method}
\alias{nullModel,missing,numeric-method}
\alias{nullModel,missing,factor-method}
\alias{nullModel,formula,missing-method}
\alias{nullModel,formula,data.frame-method}
\title{Create Null Model for Association Test}
\description{
  Method for creating a null model that can be used for
  association testing using \code{assocTest}
}
\usage{
\S4method{nullModel}{formula,data.frame}(X, y, data,
          type=c("automatic", "logistic", "linear", "bernoulli"),
          n.resampling=0,
          type.resampling=c("bootstrap", "permutation"),
          adj=c("automatic", "none", "force"), adjExact=FALSE,
          n.resampling.adj=10000, checkData=TRUE)
\S4method{nullModel}{formula,missing}(X, y, data,
          type=c("automatic", "logistic", "linear", "bernoulli"),
          n.resampling=0,
          type.resampling=c("bootstrap", "permutation"),
          adj=c("automatic", "none", "force"), adjExact=FALSE,
          n.resampling.adj=10000, checkData=TRUE)
\S4method{nullModel}{matrix,numeric}(X, y,
          type=c("automatic", "logistic", "linear"), ...)
\S4method{nullModel}{matrix,factor}(X, y,
          type=c("automatic", "logistic", "linear"), ...)
\S4method{nullModel}{missing,numeric}(X, y,
          type=c("automatic", "logistic", "linear", "bernoulli"),
          ...)
\S4method{nullModel}{missing,factor}(X, y,
          type=c("automatic", "logistic", "linear", "bernoulli"),
          ...)
}
\arguments{
  \item{X}{a formula or matrix}
  \item{y}{if the formula interface is used, \code{y} can be used to
    pass a data frame with the table in which both covariates and
    traits are contained (alternatively, the \code{data} argument
    can be used for that purpose). The other methods (if \code{X} is
    not a formula) expect \code{y} to be the trait vector.
    Trait vectors can either be numeric vectors or a factor with
    two levels (see details below).}
  \item{data}{for consistency with standard R methods from the
    \pkg{stats} package, the data frame can also be passed to
    \code{nullModel} via the \code{data} argument. In this case, the
    \code{y} must be empty. If \code{y} is specified, \code{data} is
    ignored.}
  \item{type}{type of model to train (see details below)}
  \item{n.resampling}{number of null model residuals to sample; set to
    zero (default) to turn resampling off; resampling is not supported
    for plain trait vectors without covariates}
  \item{type.resampling}{method how to sample null model residuals;
    the choice \dQuote{permutation} refers to simple random permutations
    of the model's residuals. If \dQuote{bootstrap} is chosen (default),
    the following strategy is applied for linear models (continuous
    trait): residuals are sampled as normally distributed values
    with mean 0 and the same standard deviation as the model's
    residuals. For logistic models (binary trait), the choice
    \dQuote{bootstrap} selects the same bootstrapping method that is
    implemented in the \pkg{SKAT} package.}
  \item{adj}{whether or not to use small sample correction for
    logistic models (binary trait with covariates). The choice
    \dQuote{none} turns off small sample correction. If \dQuote{force}
    is chosen, small sample correction is turned on unconditionally.
    If \dQuote{automatic} is chosen (default), small sample correction
    is turned on if the number of samples does not exceed 2,000. This
    argument is ignored for any type of model except \dQuote{logistic}
    and small sample correction is switched off.}
  \item{adjExact}{in case small sample correction is switched on (see
    above), this argument indicates whether or not the exact
    square root of the matrix \eqn{P_0} should be pre-computed (see
    Subsection 9.5 of the package vignette). The default is
    \code{FALSE}. This argument
    is ignored if small sample correction is not switched on.}
  \item{n.resampling.adj}{number of null model residuals to sample for
    the adjustment of higher moments; ignored if small sample correction
    is switched off.}
  \item{checkData}{if \code{FALSE}, only a very limited set of input
    checks is performed. The purpose of this option is to save
    computational effort for repeated input checks if the function is
    called from a function that has already performed input checks.
    The default is \code{TRUE}.
    Only change to \code{FALSE} if you know what you are doing!}
  \item{...}{all other parameters are passed on to the \code{nullModel}
    method with signature \code{formula,data.frame}.}
}
\details{
  The \pkg{podkat} package assumes a mixed model in which the
  trait under investigation depends both on covariates (if
  any) and the genotype. The \code{nullModel} method models the relationship
  between the trait and the covariates (if any) without taking the genotype
  into account, which corresponds to the null assumption that the
  trait and the genotype are independent. Therefore, we speak of
  \emph{null models}.
  The following types of models are presently
  available:
  \describe{
    \item{Linear model (type \dQuote{linear}):}{a linear
      model is trained for a continuous trait and a
      given set of covariates (if any); this is done by standard
      linear regression using the \code{\link[stats]{lm}} function.}
    \item{Logistic linear model (type \dQuote{logistic}):}{a generalized linear
      model is trained for a binary trait and a
      given set of covariates (if any); this is done by
      logistic regression using the \code{\link[stats]{glm}} function.}
    \item{Bernoulli-distributed trait (type \dQuote{bernoulli}):}{a binary
      trait without covariates is interpreted as instances of
      a simple Bernoulli process with \eqn{p} being the relative
      frequencies 1's/cases.}
  }
  The \code{type} argument can be used to select the type of model,
  where the following restrictions apply:
  \itemize{
    \item{For linear models, the trait vector must be
      numerical. Factors/factor columns are not accepted.}
    \item{For logistic models and Bernoulli-distributed traits, both numerical
      vectors and factors are acceptable. In any case, only 0's (controls)
      and 1's (cases) are accepted. Furthermore, \code{nullModel} quits
      with an error if the trait shows no variation. In other words,
      trait vectors that only contain 0's or only contain 1's are not
      accepted (as association testings makes little sense for such
      traits anyway).}
  }
  The following interfaces are available to specify the traits and the
  covariates (if any):
  \describe{
    \item{Formula interface:}{the first argument \code{X} 
      can be a formula that specifies the trait vector/column, the
      covariate matrix/columns (if any), and the intercept (if any).
      If neither the \code{y} argument nor the \code{data} argument
      is specified, \code{nullModel} searches the environment from
      which the function has been called. This interface is largely
      analogous to the functions \code{\link[stats]{lm}} and
      \code{\link[stats]{glm}}.}
    \item{Trait vector without covariates:}{if the \code{X} argument
      is omitted and \code{y} is a numeric vector or factor, \code{y} is
      interpreted as trait vector, and a null model is created from
      \code{y} without covariates. Linear and logistic models are
      trained with an intercept. For type \dQuote{bernoulli}, the
      trait vector is written to the output object as is.}
    \item{Trait vector plus covariate matrix:}{if the \code{X} argument
      is a matrix and \code{y} is a numeric vector or factor, \code{y} is
      interpreted as trait vector and \code{X} is interpreted as
      covariate matrix. In this case, linear and logistic models are
      trained as (generalized) linear regressors that predict the
      trait from the covariates plus an intercept. The type
      \dQuote{bernoulli} is not available for this variant, since
      this type of model cannot consider covariates.}
  }
  
  All \code{nullModel} methods also support the choice
  \code{type="automatic"}. In this case, \code{nullModel} guesses the
  most reasonable type of model in the following way: If the trait
  vector/column is a factor or a numeric vector containing only 0's
  and 1's (where both values must be present, as
  noted above already), the trait is supposed to be binary and
  the type \dQuote{logistic} is assumed, unless the following conditions
  are satisfied:
  \enumerate{
    \item{The number of samples does not exceed 100.}
    \item{No intercept and no covariates
      have been specified. This condition can be met by supplying an empty model
      to the formula interface (e.g. \code{y ~ 0}) or by supplying the
      trait vector as argument \code{y} while omitting \code{X}.}
  }
  If these two conditions are fulfilled for a binary trait,
  \code{nullModel} chooses the type \dQuote{bernoulli}. If the trait is
  not binary and the trait vector/column is numeric,
  \code{nullModel} assumes type \dQuote{linear}.

  For consistency with the \pkg{SKAT} package, the \pkg{podkat} package
  also offers \emph{resampling}, i.e. a certain number of vectors of
  residuals are sampled according to the null model. This can be done
  when training the null model by setting the \code{n.resampling}
  parameter (number of residual vectors that are sampled) to a value
  larger than 0. Then, when association testing is performed, p-values
  are computed also for all these sampled residuals, and an additional
  estimated p-value is computed as the relative frequency of
  p-values of sampled residuals that are at least as significant as
  the test's p-value. The procedure to sample residuals is controlled
  with the \code{type.resampling} argument (see above).
  
  For logistic models (type \dQuote{logistic}), \code{\link{assocTest}}
  offers the small sample correction as introduced by
  \cite{Lee et al. (2012)}. If the adjustment of higher moments should
  be applied, some preparations need to be
  made already when training the null model. Which preparations are
  carried out, can be controlled by the arguments \code{adj},
  \code{adjExact}, \code{n.resampling.adj}, and \code{type.resampling}
  (see descriptions of arguments above and Subsection 9.5 of the package
  vignette).

  If any missing values are found in the trait vector/column or the
  covariate matrix/columns, the respective samples are omitted from the
  resulting model (which is the standard behavior of
  \code{\link[stats]{lm}} and \code{\link[stats]{glm}} anyway).
  The indices of the omitted samples are stored in the \code{na.omit}
  slot of the returned \code{\linkS4class{NullModel}} object.
}
\value{
  returns a \code{\linkS4class{NullModel}} object
}
\author{Ulrich Bodenhofer}
\references{\url{https://github.com/UBod/podkat}

Lee, S., Emond, M. J., Bamshad, M. J., Barnes, K. C., Rieder, M. J.,
Nickerson, D. A., NHLBI Exome Sequencing Project - ESP Lung Project
Team, Christiani, D. C., Wurfel, M. M., and Lin, X. (2012)
Optimal unified approach for rare-variant association testing with
application to small-sample case-control whole-exome sequencing studies.
\emph{Am. J. Hum. Genet.} \bold{91}, 224-237. DOI:
\doi{10.1016/j.ajhg.2012.06.007}.
}
\seealso{\code{\linkS4class{NullModel}}, \code{\link[stats]{lm}},
  \code{\link[stats]{glm}}}
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

## alternatively, use the interface that only supplies the
## trait vector
model <- nullModel(y=pheno[1:100, ]$y)
model
}
\keyword{methods}
