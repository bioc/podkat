computePvalues <- function(W, x, method, symmetric=FALSE, acc=1.e-6)
{
    if (symmetric)
        lambda <- eigen(W, symmetric=TRUE, only.values=TRUE)$values
    else
        lambda <- Re(eigen(W, only.values=TRUE)$values)

    lambda <- lambda[which(lambda >= 0)]
    lambda <- lambda[which(lambda > mean(lambda) / 1.e6)]

    if (length(lambda) == 0)
    {
        warning("no positive lambdas available => setting p-values to 1",
                call.=FALSE)
        return(rep(1, length(x)))
    }

    if (method == "davies")
    {
        p.values <- .Call("davies", lambda, x, acc, PACKAGE="podkat")

        neg <- which(p.values <= 0 | p.values > 1 | !is.numeric(p.values) |
                     is.nan(p.values))

        if (length(neg) > 0)
        {
            p.values[neg] <- .Call("liuMod", lambda, x[neg], PACKAGE="podkat")
            meth <- rep("davies", length(x))
            meth[neg] <- "liu.mod"
            attr(p.values, "method") <- meth
        }
    }
    else if (method == "liu")
        p.values <- .Call("liu", lambda, x, PACKAGE="podkat")
    else
        p.values <- .Call("liuMod", lambda, x, PACKAGE="podkat")

    p.values
}

computePvaluesAdj <- function(W, kernel, x, prob, Q.sim=numeric(), method)
{
    if (kernel)
    {
        e <- eigen(W, symmetric=TRUE)

        lambda <- pmax(e$values, 0)
        U <- e$vectors

        sel <- which(lambda > mean(lambda) / 1.e6)
        lambda <- lambda[sel]
        U <- U[, sel]

        if (length(sel) == 0)
        {
            warning("no positive lambdas => resorting to bootstrap estimates",
                    call.=FALSE)
            lambda <- numeric()
        }
    }
    else
    {
        if (ncol(W) == 1)
        {
            lambda <- sum(W^2)
            U <- W / sqrt(lambda)
        }
        else
        {
            res <- try(svd(W), silent=TRUE)

            if (class(res) == "try-error")
            {
                warning("SVD error => resorting to bootstrap estimates",
                        call.=FALSE)
                lambda <- numeric()
            }
            else
            {
                lambda <- res$d^2
                U <- res$u

                sel <- which(lambda > mean(lambda) / 1.e6)
                lambda <- lambda[sel]
                U <- U[, sel]

                if (length(sel) == 0)
                {
                    warning("no positive lambdas => ",
                            "resorting to bootstrap estimates", call.=FALSE)
                    lambda <- numeric()
                }
            }
        }
    }

    if (length(lambda) > 0)
    {
        muQ <- sum(lambda)

        fac <- (3 * prob^2 - 3 * prob + 1) / (prob * (1 - prob)) - 3
        U2 <- U^2

        C <- crossprod(U2, U2 * fac)
        diag(C) <- diag(C) + 2

        varQ <- drop(crossprod(lambda, C) %*% lambda)

        if (length(Q.sim) > 0)
        {
            if (method == "unbiased")
                df <- computeChiSquareParamsUnbiased(Q.sim, muQ)$df
            else if (method == "population")
                df <- computeChiSquareParamsPopulation(Q.sim)$df
            else if (method == "sample")
                df <- computeChiSquareParamsSample(Q.sim)$df
            else
                df <- computeChiSquareParamsSKAT(Q.sim)$df

            correction <- c(exact=TRUE, resampling=TRUE)
        }
        else
        {
            lambdaS <- lambda * diag(C) / sqrt(2)
            df <- sum(lambdaS^2)^2 / sum(lambdaS^4)

            correction <- c(exact=TRUE, resampling=FALSE)
       }
    }
    else if (length(Q.sim) > 0)
    {
        if (method == "population")
            re <- computeChiSquareParamsPopulation(Q.sim)
        if (method == "sample")
            re <- computeChiSquareParamsSample(Q.sim)
        else if (method == "SKAT")
            re <- computeChiSquareParamsSKAT(Q.sim)
        else
        {
            warning("no exact mean available => ",
                    "resorting to population statistics", call.=FALSE)
            re <- computeChiSquareParamsPopulation(Q.sim)
        }

        muQ  <- re$mu
        varQ <- re$var
        df   <- re$df

        correction <- c(exact=FALSE, resampling=TRUE)
    }
    else
    {
        warning("neither positive lambdas nor bootstrap estimates available ",
                "=> setting p-values to 1", call.=FALSE)
        return(rep(1, length(x)))
    }

    pValues <- pchisq((x - muQ) * sqrt(2 * df / varQ) + df, ncp=0, df=df,
                      lower.tail=FALSE)

    attr(pValues, "correction") <- correction

    pValues
}

computeChiSquareParamsUnbiased <- function(x, muQ)
{
    mu4e <-  mean((x - muQ)^4)
    varQe <- mean((x - muQ)^2)

    if (varQe > 0)
        gamma <- mu4e / varQe^2 - 3
    else
        gamma <- 0

    # kurtosis <= 0 -> use standard value 1.e6
    if (!is.numeric(gamma) || !is.finite(gamma) || gamma <= 0)
        df <- 1.e6
    else if (gamma > 1200) # identify df by skewness
        df <-  8 * varQe^3 / mean((x - muQ)^3)^2
    else # standard case: identify df by kurtosis
        df <- 12 / gamma

    list(mu=muQ, var=varQe, df=df)
}

computeChiSquareParamsPopulation <- function(x)
{
    n <- length(x)
    muQ <- mean(x)
    mu2e <- sum((x - muQ)^2)
    mu4e <- sum((x - muQ)^4)
    varQ <- mu2e / (n - 1)

    if (varQ > 0)
        gamma <- (((n + 1) * n * mu4e) / mu2e^2 - 3 * (n - 1)) *
            (n - 1) / ((n - 2) * (n - 3))
    else
        gamma <- 0

    # kurtosis <= 0 -> use standard value 1.e6
    if (!is.numeric(gamma) || !is.finite(gamma) || gamma <= 0)
        df <- 1.e6
    else if (gamma > 1200) # identify df by skewness
    {
        mu3e <- sum((x - muQ)^3)
        df <-  8 * ((n - 2) * (n - 2) * mu2e^3) /
            (n * n * (n - 1) * mu3e^2)
    }
    else # standard case: identify df by kurtosis
        df <- 12 / gamma

    list(mu=muQ, var=varQ, df=df)
}

computeChiSquareParamsSample <- function(x)
{
    muQ <- mean(x)
    varQe <- mean((x - muQ)^2)
    mu4e <-  mean((x - muQ)^4)

    if (varQe > 0)
        gamma <- mu4e / varQe^2 - 3
    else
        gamma <- 0

    # kurtosis <= 0 -> use standard value 1.e6
    if (!is.numeric(gamma) || !is.finite(gamma) || gamma <= 0)
        df <- 1.e6
    else if (gamma > 1200) # identify df by skewness
        df <-  8 * varQe^3 / mean((x - muQ)^3)^2
    else # standard case: identify df by kurtosis
        df <- 12 / gamma

    list(mu=muQ, var=varQe, df=df)
}

computeChiSquareParamsSKAT <- function(x)
{
    muQ <- mean(x)
    mu4e <-  mean((x - muQ)^4)
    varQe <- var(x)

    if (varQe > 0)
        gamma <- mu4e / varQe^2 - 3
    else
        gamma <- 0

    # kurtosis <= 0 -> use standard value 1.e6
    if (!is.numeric(gamma) || !is.finite(gamma) || gamma <= 0)
        df <- 1.e6
    else if (gamma > 1200) # identify df by skewness
        df <-  8 * varQe^3 / mean((x - muQ)^3)^2
    else # standard case: identify df by kurtosis
        df <- 12 / gamma

    list(mu=muQ, var=varQe, df=df)
}
