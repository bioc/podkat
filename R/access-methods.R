setMethod("[",
          signature(x="GenotypeMatrix", i="index", j="missing", drop="missing"),
          function(x, i, j, ..., drop=TRUE)
          {
              addArg <- list(...)

              recomputeMAF <- FALSE
              ploidy <- 2

              if (!is.null(addArg$recomputeMAF))
              {
                  recomputeMAF <- addArg$recomputeMAF

                  recomputeMAF <- check.arg.logical(recomputeMAF)

                  if (!is.null(addArg$ploidy))
                      ploidy <- check.arg.integer(ploidy, 1, FALSE, NA, FALSE)
              }

              if (recomputeMAF)
                  xn <- genotypeMatrix(Z=x[i, , drop=FALSE], pos=x@variantInfo,
                                       ploidy=ploidy, MAF.action="ignore")
              else
              {
                  xn <- as(x[i, , drop=FALSE], "GenotypeMatrix")
                  xn@variantInfo <- x@variantInfo
              }

              xn
          })

setMethod("[",
          signature(x="GenotypeMatrix", i="missing", j="index", drop="missing"),
          function(x, i, j, ..., drop=TRUE)
          {
              xn <- as(x[, j, drop=FALSE], "GenotypeMatrix")
              xn@variantInfo <- x@variantInfo[j]

              xn
          })

setMethod("[",
          signature(x="GenotypeMatrix", i="index", j="index", drop="missing"),
          function(x, i, j, ..., drop=TRUE)
          {
              addArg <- list(...)

              recomputeMAF <- FALSE
              ploidy <- 2

              if (!is.null(addArg$recomputeMAF))
              {
                  recomputeMAF <- addArg$recomputeMAF

                  recomputeMAF <- check.arg.logical(recomputeMAF)

                  if (!is.null(addArg$ploidy))
                      ploidy <- check.arg.integer(ploidy, 1, FALSE, NA, FALSE)
              }

              if (recomputeMAF)
                  xn <- genotypeMatrix(Z=x[i, j, drop=FALSE],
                                       pos=x@variantInfo[j],
                                       ploidy=ploidy, MAF.action="ignore")
              else
              {
                  xn <- as(x[i, j, drop=FALSE], "GenotypeMatrix")
                  xn@variantInfo <- x@variantInfo[j]
              }

              xn
          })

setMethod("variantInfo", signature(x="GenotypeMatrix"),
          function(x)
          {
              out <- x@variantInfo

              if (length(colnames(x)) > 0)
                  names(out) <- colnames(x)

              out
          })

setMethod("variantInfo", signature(x="missing"),
          function(x) as(GRanges(), "VariantInfo"))

setMethod("variantInfo", signature(x="GRanges"),
          function(x) as(x, "VariantInfo"))

setMethod("MAF", signature(x="GenotypeMatrix"),
          function(x) mcols(x@variantInfo)$MAF)

setMethod("MAF", signature(x="VariantInfo"),
          function(x) mcols(x)$MAF)

setMethod("[",
          signature(x="NullModel", i="index", j="missing", drop="missing"),
          function(x, i, j, ..., drop=TRUE)
          {
              if (length(i) != length(x) ||
                  length(unique(i)) != length(i))
                  stop("'i' is not a valid permutation of samples", call.=FALSE)

              x@na.omit <- integer()
              x@residuals <- x@residuals[i]

              if (all(dim(x@model.matrix) > 0))
                  x@model.matrix <- x@model.matrix[i, , drop=FALSE]

              if (all(dim(x@res.resampling) >= 1))
                  x@res.resampling <- x@res.resampling[i, , drop=FALSE]

              if (length(x@variance) > 1)
                  x@variance <- x@variance[i]
              else
                  x@variance <- x@variance

              if (length(x@prob) > 1)
                  x@prob <- x@prob[i]
              else
                  x@prob <- x@prob

              if (all(dim(x@res.resampling.adj) >= 1))
                  x@res.resampling.adj <- x@res.resampling.adj[i, , drop=FALSE]

              if (ncol(x@P0sqrt) > 0)
                  x@P0sqrt <- x@P0sqrt[i, i]

              x
          })


setMethod("residuals", signature(object="NullModel"),
          function(object, ...) object@residuals)


setMethod("names", signature(x="NullModel"),
          function(x) names(x@residuals))


setMethod("coefficients", signature(object="NullModel"),
          function(object, ...) object@coefficients)


setMethod("length", signature(x="NullModel"), function(x) length(x@residuals))
