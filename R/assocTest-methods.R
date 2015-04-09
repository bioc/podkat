assocTest.dgCMatrix <- function(Z, model, kernel, pos, weights, width,
                                method, adj, pValueLimit)
{
    ## forward computations to specialized functions
    if (substr(kernel, 1, 6) == "linear" && model@type != "bernoulli")
        res <- assocTest.linkernel(Z=Z, model=model, kernel=kernel,
                                   weights=weights, pos=pos, width=width,
                                   method=method, adj=adj)
    else
    {
        ## compute kernel matrix
        K <- computeKernel(Z, kernel, weights=weights, pos=pos, width=width)

        ## run assocation test on kernel matrix
        res <- assocTest.kernelMatrix(Z=K, model=model,
                                      method=method, adj=adj,
                                      pValueLimit=pValueLimit)

        ## reconstruct dim info
        res@dim <- dim(Z)
    }

    ## reconstruct kernel info
    res@kernel <- kernel
    res@width <- width

    if (!is.null(weights))
        res@weights <- weights

    if (length(rownames(Z)) > 0)
        res@samples <- rownames(Z)

    res
}


assocTest.GenotypeMatrix <- function(Z, model, ranges,
                                     kernel=c("linear.podkat",
                                              "localsim.podkat",
                                              "quadratic.podkat",
                                              "linear.SKAT",
                                              "localsim.SKAT",
                                              "quadratic.SKAT"),
                                     width=1000, weights=NULL,
                                     weightFunc=betaWeights(),
                                     method=NULL,
                                     adj=c("automatic", "none", "force"),
                                     pValueLimit=0.05)
{
    ## input checks
    kernel <- match.arg(kernel)
    adj    <- match.arg(adj)
    pValueLimit <- check.arg.numeric(pValueLimit, 0, TRUE, 1, FALSE)

    if (model@type != "logistic" ||
        (adj == "automatic" && length(model) > 2000))
        adj <- "none"

    if (adj == "none")
        method <- match.arg(method, c("davies", "liu.mod", "liu"))
    else
        method <- match.arg(method, c("unbiased", "population", "sample",
                                      "SKAT"))

    if (any(dim(Z) == 0)) ## empty matrix => return default result
    {
        warning("no data in genotype matrix 'Z'", call.=FALSE)

        if (missing(ranges))
        {
            res <- new("AssocTestResult",
                       type=model@type,
                       kernel=kernel,
                       dim=dim(Z),
                       weights=numeric(),
                       width=width,
                       method="",
                       Q=0,
                       p.value=1,
                       Q.resampling=numeric(),
                       p.value.resampling=numeric(),
                       p.value.resampled=1,
                       call=deparse(sys.call(-1)))

            if (ncol(model@res.resampling) > 0)
            {
                res@Q.resampling <- rep(0, ncol(model@res.resampling))
                res@p.value.resampling <- rep(1, ncol(model@res.resampling))
            }
        }
        else
        {
            res <- as(ranges, "AssocTestResultRanges")

            res@type <- model@type
            res@kernel <- kernel
            res@call <- deparse(sys.call(-1))
            res@weights=NULL
            res@width=width

            if (length(rownames(Z)) > 0)
                ranges@samples <- rownames(Z)

            if (length(res) > 0)
            {
                mcols(res)$n <- 0
                mcols(res)$Q <- 0
                mcols(res)$p.value <- 1
            }
        }

        return(res)
    }

    if (is.null(weights))
    {
        ## compute weights
        if (!is.null(weightFunc))
        {
            if (!is.function(weightFunc) || length(formals(weightFunc)) < 1)
                stop("'weightFunc' must be NULL or function with at least one ",
                     "argument", call.=FALSE)

            weights <- try(weightFunc(mcols(Z@variantInfo)$MAF))

            if (class(weights) == "try-error" || !is.numeric(weights) ||
                length(weights) != ncol(Z) || any(weights < 0))
                stop("invalid weighting function 'weightFunc'", call.=FALSE)
        }
    }
    else if (is.numeric(weights))
    {
        if (length(weights) != ncol(Z))
            stop("weight vector must be NULL or exactly as long as ",
                 "number of variants", call.=FALSE)
    }
    else
        stop("invalid 'weights' argument", call.=FALSE)

    if (!is.numeric(width) || length(width) != 1 || width <= 0)
        stop("invalid 'width' argument", call.=FALSE)

    ## check for correct order of samples
    if (length(rownames(Z)) > 0 && length(names(model@residuals)) > 0)
    {
        sel <- match(names(model@residuals), rownames(Z))

        if (any(is.na(sel)))
            stop("not all samples of null model contained in genotype matrix",
                 call.=FALSE)

        Z <- Z[sel, ]
    } ## if no names are present, assume that order is correct
    else if (nrow(Z) == (length(model) + length(model@na.omit)))
    {
        if (length(model@na.omit) > 0)
            Z <- Z[-model@na.omit, ]
    }
    else if (nrow(Z) != length(model))
        stop("dimension of genotype matrix 'Z' does not fit to null model",
             call.=FALSE)

    ## variant 1: no 'ranges' argument => single test
    if (missing(ranges))
    {
       ## only regions from the same chromosome are allowed (otherwise
        ##     position weighting is messed up)
        if (length(unique(seqnames(Z@variantInfo))) > 1)
            stop("genotype matrix 'Z' contains variants from more than one ",
                 "chromosome", call.=FALSE)

        ## run test
        res <- assocTest.dgCMatrix(Z=as(Z, "dgCMatrix"), model=model,
                                   kernel=kernel, pos=start(Z@variantInfo),
                                   weights=weights, width=width,
                                   method=method, adj=adj,
                                   pValueLimit=pValueLimit)

        res@call <- deparse(sys.call(-1))

        res
    }
    else ## multiple windows along 'ranges' argument
    {
        n <- length(ranges)

        chr <- seqnames(Z@variantInfo)
        pos <- start(Z@variantInfo)

        out <- data.frame(n=rep(0, n), Q=rep(0, n), p.value=rep(1, n),
                          p.value.resampled=rep(1, n))

        env <- environment()

        ## function running individual tests
        assocTest.singleWindow <- function(i)
        {
            ## select variants in selected region
            subsel <- which(as.logical(chr ==
                                       as.character(seqnames(ranges)[i]) &
                                       pos >= start(ranges)[i] &
                                       pos <= end(ranges)[i]))

            if (length(subsel) > 0)
            {
                res <- assocTest.dgCMatrix(Z=Z[, subsel],
                                           model=model, kernel=kernel,
                                           pos=pos[subsel],
                                           weights=weights[subsel], width=width,
                                           method=method, adj=adj,
                                           pValueLimit=pValueLimit)

                env$out[i, ] <- c(res@dim[2], res@Q, res@p.value,
                                  res@p.value.resampled)
            }

            invisible(NULL)
        }

        ## run individual tests
        lapply(seq(1, length.out=n), assocTest.singleWindow)

        mcols(ranges) <- out

        ranges <- as(ranges, "AssocTestResultRanges", strict=FALSE)

        ranges@type    <- model@type
        ranges@kernel  <- kernel
        ranges@weights <- weights
        ranges@width   <- width
        ranges@call    <- deparse(sys.call(-1))

        if (length(rownames(Z)) > 0)
            ranges@samples <- rownames(Z)

        if (ncol(model@res.resampling) == 0)
            mcols(ranges)$p.value.resampled <- NULL

        ranges
    }
}

setMethod("assocTest",
          signature(Z="GenotypeMatrix", model="NullModel"),
          assocTest.GenotypeMatrix)


assocTest.linkernel <- function(Z, model, kernel, weights, pos, width=1000,
                                method, adj)
{
    ## scale matrix by weights
    if (!is.null(weights))
        Z <- as(Z, "dgCMatrix") %*% Diagonal(length(weights), weights)

    kernel <- substr(kernel, 8, nchar(kernel))

    ## convolve matrix with position kernel
    if (kernel == "podkat")
        Z <- as(Z, "dgCMatrix") %*%
             .Call("posKernel", as.integer(pos), as.double(width),
                   PACKAGE="podkat")

    ## compute test statistic
    Qres <- crossprod(model@residuals, Z)
    Q <- drop(tcrossprod(Qres)) / 2

    ## compute resampling test statistics
    if (ncol(model@res.resampling) > 0)
    {
        Qresres <- crossprod(model@res.resampling, Z)
        Q.resampling <- rowSums(Qresres^2) / 2
    }
    else
        Q.resampling <- numeric()

    if (model@type == "linear") ## linear model
    {
        Q <- Q / model@variance
        Q.resampling <- Q.resampling / model@variance

        if (nrow(Z) > ncol(Z)) ## choose faster variant
        {
            Zpi <- crossprod(model@model.matrix, Z)
            W <- crossprod(Z) - crossprod(Zpi, model@inv.matrix) %*% Zpi
        }
        else
        {
            P0Z <- Z - (model@model.matrix %*% model@inv.matrix) %*%
                       crossprod(model@model.matrix, Z)
            W <- tcrossprod(P0Z)
        }

        ## compute p-values
        p.values <- computePvalues(W / 2, x=c(Q, Q.resampling),
                                   method=method, symmetric=TRUE, acc=1.e-6)
    }
    else ## logistic model
    {
        if (adj != "none") ## small sample correction
        {
            ## compute resampling statistics
            if (ncol(model@res.resampling.adj) > 0)
            {
                Qsimres <- t(model@res.resampling.adj) %*% Z
                Q.sim <- rowSums(Qsimres^2) / 2
            }
            else
                Q.sim <- numeric()

            if (nrow(model@P0sqrt) >= 1) ## exact variant
                Z1 <- model@P0sqrt %*% Z / sqrt(2)
            else ## approximate variant (like in SKAT package)
            {
                varSq <- sqrt(model@variance / 2)
                Z1 <- Z * varSq -
                    (model@model.matrix * varSq) %*%
                        model@inv.matrix %*%
                            crossprod(model@model.matrix,
                                      Z * model@variance)
            }

            ## compute p-values
            p.values <- computePvaluesAdj(Z1, kernel=FALSE,
                                          x=c(Q, Q.resampling),
                                          prob=model@prob, Q.sim=Q.sim,
                                          method=method)
        }
        else
        {
            if (nrow(Z) > ncol(Z)) ## choose faster variant
            {
                Zpi <- Z * model@variance
                XZpi <- crossprod(model@model.matrix, Zpi)

                W <- crossprod(Z, Zpi) -
                     crossprod(XZpi, model@inv.matrix) %*% XZpi

                ## compute p-values
                p.values <- computePvalues(W / 2, x=c(Q, Q.resampling),
                                           method=method, symmetric=TRUE,
                                           acc=1.e-6)
            }
            else
            {
                VX <- model@model.matrix * model@variance
                K <- tcrossprod(Z)
                W <- K * model@variance -
                     VX %*% model@inv.matrix %*%
                     crossprod(VX, K)

                ## compute p-values
                p.values <- computePvalues(W / 2, x=c(Q, Q.resampling),
                                           method=method, symmetric=FALSE,
                                           acc=1.e-6)
            }
        }
    }

    if (is.null(p.values))
        stop("error determining null distribution", call.=FALSE)

    ## compute resampling p-values
    if (ncol(model@res.resampling) > 0)
        p.value.resampled <- length(which(p.values <= p.values[1])) /
                             length(p.values)
    else
        p.value.resampled <- 1

    ## put together result object
    res <- new("AssocTestResult",
               type=model@type,
               kernel=kernel,
               dim=dim(Z),
               Q=Q,
               p.value=p.values[1],
               Q.resampling=Q.resampling,
               p.value.resampling=p.values[-1],
               p.value.resampled=p.value.resampled)

    if (!is.null(attr(p.values, "correction")))
        res@correction <- attr(p.values, "correction")

    if (length(attr(p.values, "method")) > 0)
    {
        meth <- attr(p.values, "method")
        res@method <- list(Q=meth[1], Q.resampling=meth[-1])
    }
    else
        res@method <- method

    res
}


assocTest.kernelMatrix <- function(Z, model, method=NULL,
                                   adj=c("automatic", "none", "force"),
                                   pValueLimit=0.05)
{
    ## input checks
    adj <- match.arg(adj)
    pValueLimit <- check.arg.numeric(pValueLimit, 0, TRUE, 1, FALSE)

    if (model@type != "logistic" ||
        (adj == "automatic" && length(model) > 2000))
        adj <- "none"

    if (adj == "none")
        method <- match.arg(method, c("davies", "liu.mod", "liu"))
    else
        method <- match.arg(method, c("unbiased", "population", "sample",
                                      "SKAT"))

    if (ncol(Z) != nrow(Z))
        stop("kernel matrix 'Z' must be quadratic", call.=FALSE)

    ## check for correct order of samples
    if (length(rownames(Z)) > 0 && length(names(model@residuals)) > 0)
    {
        sel <- match(names(model@residuals), rownames(Z))

        if (any(is.na(sel)))
            stop("not all samples of null model contained in genotype matrix",
                 call.=FALSE)

        Z <- Z[sel, sel]
    } ## if no names are present, assume that order is correct
    else if (nrow(Z) == (length(model) + length(model@na.omit)))
    {
        if (length(model@na.omit) > 0)
            Z <- Z[-model@na.omit, -model@na.omit]
    }
    else if (nrow(Z) != length(model))
        stop("dimension of genotype matrix 'Z' does not fit to null model",
             call.=FALSE)


    else if (nrow(Z) != length(model))
        stop("numbers of samples do not match between kernel matrix and ",
             "null model", call.=FALSE)

    if (!is.numeric(pValueLimit) || length(pValueLimit) != 1 ||
        pValueLimit <= 0 || pValueLimit > 1)
        stop("'pValueLimit' must be numeric value above 0 and not larger ",
             "than 1", call.=FALSE)

    n <- ncol(Z)

    Q <- drop(model@residuals %*% Z %*% model@residuals) / 2

    if (ncol(model@res.resampling) > 0)
        Q.resampling <- colSums((Z %*% model@res.resampling) *
                                model@res.resampling) / 2
    else
        Q.resampling <- numeric()

    if (model@type == "linear") ## linear model
    {
        Q <- Q / model@variance
        Q.resampling <- Q.resampling / model@variance
        P0K <- Z - (model@model.matrix %*% model@inv.matrix) %*%
                    crossprod(model@model.matrix, Z)
        W <- P0K - (P0K %*% model@model.matrix) %*%
                   tcrossprod(model@inv.matrix, model@model.matrix)

        ## compute p-values
        p.values <- computePvalues(W / 2, x=c(Q, Q.resampling),
                                   method=method, symmetric=TRUE, acc=1.e-6)
    }
    else if (model@type == "logistic") ## logistic model
    {
        if (adj != "none") ## small sample correction
        {
            if (ncol(model@res.resampling.adj) > 0)
                 Q.sim <- colSums((Z %*% model@res.resampling.adj) *
                                 model@res.resampling.adj) / 2
            else
                Q.sim <- numeric()

            if (nrow(model@P0sqrt) >= 1) ## exact adjustment
                W <- model@P0sqrt %*% Z %*% model@P0sqrt
            else ## approximate variant (like in SKAT package)
            {
                varSq <- sqrt(model@variance)
                VX <- model@model.matrix * model@variance
                W <- (Z * varSq) - (model@model.matrix * varSq) %*%
                      model@inv.matrix %*% crossprod(VX, Z)
                W <- sweep(W, MARGIN=2, STATS=varSq, FUN="*") -
                     (W %*% VX) %*% tcrossprod(model@inv.matrix,
                                               model@model.matrix * varSq)
            }

            ## compute p-values
            p.values <- computePvaluesAdj(W / 2, kernel=TRUE,
                                          x=c(Q, Q.resampling),
                                          prob=model@prob, Q.sim=Q.sim,
                                          method=method)
        }
        else
        {
            VX <- model@model.matrix * model@variance
            W <- Z * model@variance -
                 VX %*% model@inv.matrix %*% crossprod(VX, Z)

            ## compute p-values
           p.values <- computePvalues(W / 2, x=c(Q, Q.resampling),
                                      method=method, symmetric=FALSE, acc=1.e-6)
        }
    }
    else ## exact mixture-of-Bernoulli test
    {
        p.values <- .Call("computeExactBernoulliPvalue", as.double(Q),
                          Z, as.double(model@prob), as.double(pValueLimit),
                          PACKAGE="podkat")

        if (p.values >= pValueLimit)
            p.values <- 1
    }

    if (is.null(p.values))
        stop("error determining null distribution", call.=FALSE)

    ## compute resampling p-values
    if (ncol(model@res.resampling) > 0)
        p.value.resampled <- length(which(p.values <= p.values[1])) /
                             length(p.values)
    else
        p.value.resampled <- 1

    ## put together result object
    res <- new("AssocTestResult",
               type=model@type,
               kernel="user-defined",
               Q=Q,
               p.value=p.values[1],
               Q.resampling=Q.resampling,
               p.value.resampling=p.values[-1],
               p.value.resampled=p.value.resampled)

    if (length(rownames(Z)) > 0)
        res@samples <- rownames(Z)

    if (length(attr(p.values, "method")) > 0)
    {
        meth <- attr(p.values, "method")
        res@method <- list(Q=meth[1], Q.resampling=meth[-1])
    }
    else
        res@method <- method

    res@call <- deparse(sys.call(-1))

    res
}

setMethod("assocTest",
          signature(Z="matrix", model="NullModel"),
          assocTest.kernelMatrix)


## worker routine (does not do any input checks)

assocTest.TabixFileWorker <- function(ranges, file, model, kernel, subset,
                                      sex, width, noIndels, onlyPass, na.limit,
                                      MAF.limit, na.action, MAF.action,
                                      weightFunc, method, adj,
                                      pValueLimit=pValueLimit)
{
    n <- length(ranges)

    out <- data.frame(n=rep(0, n), Q=rep(0, n), p.value=rep(1, n),
                      p.value.resampled=rep(1, n))

    if (n == 0)
        return(out)

    rranges <- reduce(ranges) ## read whole batch at once

    ## read VCF file
    vcf <- .vcfScan(file, seqnames=as.character(seqnames(rranges)),
                    start=start(rranges), end=end(rranges), subset=subset,
                    noIndels=noIndels, onlyPass=onlyPass, na.limit=na.limit,
                    MAF.limit=MAF.limit, na.action=na.action,
                    MAF.action=MAF.action, sex=sex)

    if (is.null(vcf)) ## no variants => exit with empty result
        return(out)

    if (is.character(model))
    {
        if (exists("model", envir=globalenv()))
            model <- .GlobalEnv$model
        else
            stop("'model' not available on client node", call.=FALSE)
    }

    chr <- vcf$seqnames
    pos <- vcf$pos
    Z <- vcf$GT

    ## compute weight vector
    if (!is.null(weightFunc))
    {
        weights <- try(weightFunc(vcf$MAF))

        if (class(weights) == "try-error" || !is.numeric(weights) ||
            length(weights) != length(vcf$MAF) || any(weights < 0))
        {
            warning("error computing weight vector => using unweighted test",
                    call.=FALSE)

            weights <- NULL
        }
    }
    else
        weights <- NULL

    env <- environment()

    ## function running individual tests
    assocTest.singleWindow <- function(i)
    {
        ## select variants in selected region
        subsel <- which(as.logical(chr == as.character(seqnames(ranges)[i]) &
                                   pos >= start(ranges)[i] &
                                   pos <= end(ranges)[i]))

        if (length(subsel) > 0)
        {
            res <- assocTest.dgCMatrix(Z=Z[, subsel, drop=FALSE], model=model,
                                       kernel=kernel, pos=pos[subsel],
                                       weights=if (is.null(weights)) NULL
                                               else weights[subsel],
                                       width=width, method=method, adj=adj,
                                       pValueLimit=pValueLimit)

            env$out[i, ] <- c(res@dim[2], res@Q, res@p.value,
                              res@p.value.resampled)
        }

        invisible(NULL)
    }

    ## run individual tests
    lapply(1:length(ranges), assocTest.singleWindow)

    out
}


assocTest.TabixFile <- function(Z, model, ranges,
                                kernel=c("linear.podkat", "localsim.podkat",
                                         "quadratic.podkat", "linear.SKAT",
                                         "localsim.SKAT", "quadratic.SKAT"),
                                cl=NULL, nnodes=1,
                                batchSize=20, noIndels=TRUE, onlyPass=TRUE,
                                na.limit=1, MAF.limit=1,
                                na.action=c("impute.major", "omit"),
                                MAF.action=c("invert", "omit", "ignore"),
                                sex=NULL, weightFunc=betaWeights(), width=1000,
                                method=NULL,
                                adj=c("automatic", "none", "force"),
                                pValueLimit=(0.1 / length(ranges)),
                                tmpdir=tempdir(),
                                displayProgress=TRUE)
{
    ## input checks
    kernel      <- match.arg(kernel)
    adj         <- match.arg(adj)
    na.action   <- match.arg(na.action)
    MAF.action  <- match.arg(MAF.action)
    noIndels    <- check.arg.logical(noIndels)
    onlyPass    <- check.arg.logical(onlyPass)
    MAF.limit   <- check.arg.numeric(MAF.limit, 0, TRUE, 1, FALSE)
    na.limit    <- check.arg.numeric(na.limit, 0, TRUE, 1, FALSE)
    pValueLimit <- check.arg.numeric(pValueLimit, 0, TRUE, 1, FALSE)
    width       <- check.arg.numeric(width, 0, TRUE, NA, FALSE)
    nnodes      <- check.arg.integer(nnodes, 1, FALSE, NA, FALSE)

    if (!is.numeric(batchSize) || length(batchSize) < 1 ||
        length(batchSize) > 2 || any(round(batchSize) != batchSize) ||
        any(batchSize < 1) ||
        (length(batchSize) > 1 && batchSize[1] > batchSize[2]))
        stop("invalid 'batchSize' argument", call.=FALSE)

    if (length(batchSize) == 2 && batchSize[1] == batchSize[2])
        batchSize <- batchSize[1]

    if (is(ranges, "GRanges"))
    {
        if (length(ranges) == 0)
            stop("'ranges' must be a non-empty", call.=FALSE)
    }
    else if (is(ranges, "GRangesList"))
    {
        if (length(ranges) == 0)
            stop("'ranges' must be a non-empty", call.=FALSE)

        sel <- which(sapply(ranges, length) > 0)

        if (length(sel) == 0)
            stop("'ranges' does not contain non-empty components", call.=FALSE)

        ranges <- ranges[sel]
    }
    else
        stop("'ranges' must be a non-empty object of class 'GRanges' or ",
             "'GRangesList'", call.=FALSE)

    if (!is.null(cl) && !is(cl, "SOCKcluster"))
        stop("'cl' must be NULL or an object of class 'SOCKcluster'",
             call.=FALSE)

    if (model@type != "logistic" ||
        (adj == "automatic" && length(model) > 2000))
        adj <- "none"

    if (adj == "none")
    {
        method <- match.arg(method, c("davies", "liu.mod", "liu"))
        model@res.resampling.adj <- matrix(nrow=0, ncol=0)
    }
    else
        method <- match.arg(method, c("unbiased", "population", "sample",
                                      "SKAT"))

    if (!is.null(weightFunc))
    {
        if (!is.function(weightFunc) || length(formals(weightFunc)) < 1)
            stop("'weightFunc' must be NULL or function with at least one ",
                 "argument", call.=FALSE)

        res <- try(weightFunc(seq(0.0001, 0.9999, length.out=20)))

        if (class(weights) == "try-error" || !is.numeric(res) ||
            length(res) != 20 || any(res < 0))
            stop("invalid weighting function 'weightFunc'", call.=FALSE)
    }

    if (!isOpen(Z))
    {
        open(Z)
        on.exit(close(Z))
    }

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(names(sex)) == 0 && length(sex) != length(model))
            stop("length of 'sex' argument does not match number of ",
                 "samples in null model", call.=FALSE)
    }

    if (is.factor(sex))
    {
        lev <- levels(sex)

        if (length(lev) != 2 || any(lev != c("F", "M")))
            stop("invalid 'sex' argument", call.=FALSE)
    }

    ## check for data integrity and possible sub-setting
    sampleNames <- readSampleNamesFromVcfHeader(Z)

    if (length(sampleNames) == 0)
        stop("could not determine sample names from tabix file ",
             path(Z), call.=FALSE)

    ## try to determine sample names from null model
    sampleNamesNM <- names(model@residuals)

    subset <- logical()

    if (length(sampleNamesNM) > 0)
    {
        fileN <- length(sampleNames)

        subset <- (sampleNames %in% sampleNamesNM)

        if (length(which(subset)) != length(sampleNamesNM))
            stop("some samples of null model are missing in tabix file ",
                 path(Z), call.=FALSE)

        sampleNames <- sampleNames[subset]

        perm <- match(sampleNames, sampleNamesNM)

        model <- model[perm]

        if (!is.null(sex))
        {
            if (length(names(sex)) > 0)
            {
                perm <- match(sampleNames, names(sex))

                if (any(is.na(perm)))
                    stop("name mismatch between 'sex' and null model",
                         call.=FALSE)
            }

            sex <- sex[perm]
        }
    }
    else
        stop("no sample names in null model", call.=FALSE)

    if (is(ranges, "GRangesList"))
    {
        if (length(names(ranges)) > 0)
            componentNames <- names(ranges)
        else
            componentNames <- as.character(1:length(ranges))

        if (length(mcols(ranges)$doubleMales) > 0)
            sexL <- lapply(mcols(ranges)$doubleMales,
                           function(chAs) if (chAs) sex else NULL)
        else
            sexL <- lapply(1:length(ranges), function(i) NULL)

        if (displayProgress)
        {
            startTime <- Sys.time()
            message("\tstarting association test")
        }

        if (is.null(cl) && nnodes > 1)
        {
            cl <- makePSOCKcluster(nnodes)
            on.exit(stopCluster(cl))

            if (displayProgress)
                message("\tcluster with ", nnodes,
                        " nodes started (elapsed time: ",
                        format(Sys.time() - startTime, digits=3), ")")
        }

        res <- lapply(1:length(ranges),
                      function(chr)
                      {
                          out <- assocTest(Z=Z,
                                           model=model,
                                           ranges=ranges[[chr]],
                                           kernel=kernel,
                                           cl=cl,
                                           batchSize=batchSize,
                                           noIndels=noIndels,
                                           onlyPass=onlyPass,
                                           na.limit=na.limit,
                                           MAF.limit=MAF.limit,
                                           na.action=na.action,
                                           MAF.action=MAF.action,
                                           sex=sexL[[chr]],
                                           weightFunc=weightFunc,
                                           width=width,
                                           method=method,
                                           adj=adj,
                                           pValueLimit=pValueLimit)

                          if (displayProgress)
                              message("\tdone with component '",
                                      componentNames[chr], "' (elapsed time: ",
                                      format(Sys.time() - startTime, digits=3),
                                      ")")

                          out
                      })

        ranges <- do.call(c, unname(res))

        if (displayProgress)
            message("\toutput merged => finished (elapsed time: ",
                    format(Sys.time() - startTime, digits=3))
    }
    else if (length(ranges) <= batchSize[1])
        mcols(ranges) <- assocTest.TabixFileWorker(ranges=ranges, file=Z,
                                                   model=model, kernel=kernel,
                                                   subset=subset, sex=sex,
                                                   noIndels=noIndels,
                                                   onlyPass=onlyPass,
                                                   na.limit=na.limit,
                                                   MAF.limit=MAF.limit,
                                                   na.action=na.action,
                                                   MAF.action=MAF.action,
                                                   weightFunc=weightFunc,
                                                   width=width,
                                                   method=method, adj=adj,
                                                   pValueLimit=pValueLimit)
    else if (is.null(cl) && nnodes == 1)
    {
        if (length(batchSize) > 1)
        {
            warning("'batchSize' has length 2 even though 'nnodes' is 1; ",
                    "=> ignoring first value", call.=FALSE)
            batchSize <- batchSize[2]
        }

        indexList <- lapply(seq(from=1, to=length(ranges), by=batchSize),
                            function(i)
                                seq(from=i,
                                    to=min(length(ranges), i + batchSize - 1),
                                    by=1))

        func <- function(indices, ...)
            assocTest.TabixFileWorker(ranges[indices], ...)

        mcols(ranges) <- do.call(rbind,
                                 lapply(indexList, func,
                                        file=Z, model=model, kernel=kernel,
                                        subset=subset, sex=sex,
                                        noIndels=noIndels, onlyPass=onlyPass,
                                        na.limit=na.limit, MAF.limit=MAF.limit,
                                        na.action=na.action,
                                        MAF.action=MAF.action,
                                        weightFunc=weightFunc, width=width,
                                        method=method, adj=adj,
                                        pValueLimit=pValueLimit))
    }
    else
    {
        if (is.null(cl))
        {
            if (length(batchSize) == 1)
                nnodes <- min(nnodes, ceiling(length(ranges) / batchSize))

            cl <- makePSOCKcluster(nnodes)
            on.exit(stopCluster(cl))
        }
        else
            nnodes <- length(cl)

        if (length(batchSize) > 1)
            batchSize <- min(batchSize[2],
                             max(batchSize[1],
                                 floor(length(ranges) / (4 * nnodes))))

        indexList <- lapply(seq(from=1, to=length(ranges), by=batchSize),
                            function(i)
                                seq(from=i,
                                    to=min(length(ranges), i + batchSize - 1),
                                    by=1))

        ret <- clusterEvalQ(cl, library(podkat))

        func <- function(indices, ...)
            assocTest.TabixFileWorker(ranges[indices], ...)

        tmpfile <- tempfile(pattern="model", tmpdir=tmpdir, fileext=".RData")
        save(model, file=tmpfile)

        ret <- clusterCall(cl, function(x) load(file=x, envir=globalenv()),
                           tmpfile)

        mcols(ranges) <- do.call(rbind,
                                 clusterApplyLB(cl, indexList, func, file=Z,
                                                model=tmpfile, kernel=kernel,
                                                subset=subset, sex=sex,
                                                noIndels=noIndels,
                                                onlyPass=onlyPass,
                                                na.limit=na.limit,
                                                MAF.limit=MAF.limit,
                                                na.action=na.action,
                                                MAF.action=MAF.action,
                                                weightFunc=weightFunc,
                                                width=width,
                                                method=method, adj=adj,
                                                pValueLimit=pValueLimit))

        unlink(tmpfile)
    }

    ranges <- as(ranges, "AssocTestResultRanges", strict=FALSE)

    ranges@type      <- model@type
    ranges@samples   <- sampleNames
    ranges@kernel    <- kernel
    ranges@weights   <- weightFunc
    ranges@width     <- width
    ranges@vcfParams <- list(noIndels=noIndels, onlyPass=onlyPass,
                             na.limit=na.limit, MAF.limit=MAF.limit,
                             na.action=na.action, MAF.action=MAF.action)
    ranges@call      <- deparse(sys.call(-1))

    if (!is.null(sex))
        ranges@sex <- sex

    if (ncol(model@res.resampling) <= 0)
        mcols(ranges)$p.value.resampled <- NULL

    ranges
}

setMethod("assocTest",
          signature(Z="TabixFile", model="NullModel"),
          assocTest.TabixFile)


assocTest.filename <- function(Z, model, ...)
{
    res <- assocTest(TabixFile(Z), model, ...)

    res@call <- deparse(sys.call(-1))

    res
}

setMethod("assocTest",
          signature(Z="character", model="NullModel"),
          assocTest.filename)
