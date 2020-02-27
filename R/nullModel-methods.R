#
# nullModel: variant with covariate matrix and response
#
nullModel.X.y <- function(X, y,
                          type=c("automatic", "logistic", "linear"),
                          ...)
{
    type <- match.arg(type)

    if (!is.vector(y) && !is.factor(y))
        stop("response 'y' must be vector or factor", call.=FALSE)

    if (nrow(X) != length(y))
        stop("numbers of samples in 'X' and 'y' do not match", call.=FALSE)

    if (type == "automatic")
    {
        if ((is.factor(y) && length(levels(y)) == 2) ||
            length(unique(y)) == 2)
            type <- "logistic"
        else if (is.vector(y) && is.numeric(y))
            type <- "linear"
        else
            stop("response 'y' has invalid type", call.=FALSE)
    }

    if (type == "logistic")
    {
        if (is.factor(y))
        {
            yn <- names(y)
            y <- as.numeric(as.character(y))
            names(y) <- yn
        }

        if (any(!(y %in% c(0, 1))))
            stop("for 'type=\"logistic\"', response 'y' must be binary",
                 call.=FALSE)
    }
    else if (type == "linear")
    {
        if (!is.vector(y) || !is.numeric(y))
            stop("for 'type=\"linear\"', response 'y' must be numeric vector",
                 call.=FALSE)
    }
    else
        stop("invalid 'type' argument", call.=FALSE)

    out <- nullModel(y ~ X, type=type, ..., checkData=FALSE)

    out@call <- deparse(sys.call(-1))

    out
}

setMethod("nullModel", signature(X="matrix", y="numeric"), nullModel.X.y)
setMethod("nullModel", signature(X="matrix", y="factor"), nullModel.X.y)


#
# nullModel: variant only with response
#
nullModel.y <- function(X, y,
                        type=c("automatic", "logistic", "linear", "bernoulli"),
                        ...)
{
    type <- match.arg(type)

    if (type == "automatic")
    {
        if ((is.factor(y) && length(levels(y)) == 2) ||
            length(unique(y)) == 2)
        {
            if (length(y) <= 100)
                type <- "bernoulli"
            else
                type <- "logistic"
        }
        else if (is.vector(y) && is.numeric(y))
            type <- "linear"
        else
            stop("response 'y' has invalid type", call.=FALSE)
    }

    if (type == "logistic" || type == "bernoulli")
    {
        if (is.factor(y))
        {
            yn <- names(y)
            y <- as.numeric(as.character(y))
            names(y) <- yn
        }

        if (length(unique(y)) != 2 || any(!(y %in% c(0, 1))))
            stop("for 'type=\"", type, "\"', response 'y' must be binary",
                 call.=FALSE)
    }
    else if (type == "linear")
    {
        if (!is.vector(y) || !is.numeric(y))
             stop("for 'type=\"linear\"', response 'y' must be numeric vector",
                 call.=FALSE)
    }
    else
        stop("invalid 'type' argument", call.=FALSE)

    if (type == "bernoulli")
        out <- nullModel(y ~ 0, type=type, ..., checkData=FALSE)
    else
        out <- nullModel(y ~ 1, type=type, ..., checkData=FALSE)

    out@call <- deparse(sys.call(-1))

    out
}

setMethod("nullModel", signature(X="missing", y="numeric"), nullModel.y)
setMethod("nullModel", signature(X="missing", y="factor"), nullModel.y)


#
# nullModel: variant with formula
#
nullModel.formula <- function(X, y, data,
                              type=c("automatic", "logistic", "linear",
                                     "bernoulli"),
                              n.resampling=0,
                              type.resampling=c("bootstrap", "permutation"),
                              adj=c("automatic", "none", "force"),
                              adjExact=FALSE, n.resampling.adj=10000,
                              checkData=TRUE)
{
    type             <- match.arg(type)
    type.resampling  <- match.arg(type.resampling)
    adj              <- match.arg(adj)
    n.resampling     <- check.arg.integer(n.resampling, 0, FALSE, NA, FALSE)
    n.resampling.adj <- check.arg.integer(n.resampling.adj, 0, FALSE, NA, FALSE)

    if (missing(y))
    {
        y <- NULL

        if (!missing(data))
        {
            if (is.data.frame(data))
                y <- data
            else
                stop("invalid 'data' argument", call.=FALSE)
        }
    }

    if (checkData)
    {
        yr <- model.response(model.frame(X, data=y))

        if (type == "automatic")
        {
            if ((is.factor(yr) && length(levels(yr)) == 2) ||
                length(levels(factor(yr))) == 2)
            {
                type <- "logisticOrBernoulli"

                if (is.factor(yr))
                    yr <- as.numeric(as.character(yr))
            }
            else if (is.vector(yr) && is.numeric(yr))
                type <- "linear"
            else
                stop("response column has invalid type", call.=FALSE)
        }

        if (type == "logistic" || type == "bernoulli")
        {
            if (is.factor(yr))
                yr <- as.numeric(as.character(yr))

            if (length(unique(yr)) != 2 || any(!(yr %in% c(0, 1))))
                stop("for 'type=\"", type, "\"', response 'y' must be binary",
                     call.=FALSE)
        }
        else if (type == "linear")
        {
            if (!is.vector(yr) || !is.numeric(yr))
                stop("for 'type=\"linear\"', response must be binary",
                     call.=FALSE)
        }
        else if (type != "logisticOrBernoulli")
            stop("invalid 'type' argument", call.=FALSE)
    }

    mf <- model.frame(X, data=y, na.action=na.omit)

    if (length(attr(mf, "na.action")) > 0)
    {
        omitted <- as.integer(attr(mf, "na.action"))
        warning("NAs present => ", length(omitted), " samples omitted",
                call.=FALSE)
    }
    else
        omitted <- integer()

    moMat <- model.matrix(X, data=y)

    if (type == "logisticOrBernoulli")
    {
        if (ncol(moMat) == 0)
            type <- "bernoulli"
        else
            type <- "logistic"
    }

    if (type == "bernoulli")
    {
        if (ncol(moMat) > 0)
            stop("no covariates allowed for type 'exact'", call.=FALSE)
    }
    else if (ncol(moMat) == 0)
        stop("at least intercept is required for type '", type, "'",
             call.=FALSE)

    if (type == "logistic")
    {
        glmres <- glm(X, data=y, family="binomial")

        moMat <- model.matrix(glmres)

        out <- new("NullModel",
                   type=type,
                   residuals=(glmres$y - glmres$fitted.values),
                   model.matrix=moMat,
                   coefficients=glmres$coefficients,
                   na.omit=omitted,
                   n.cases=sum(glmres$y),
                   type.resampling=type.resampling,
                   prob=glmres$fitted.values,
                   variance=(glmres$fitted.values *
                             (1 - glmres$fitted.values)))

        naCoeff <- which(is.na(out@coefficients))

        if (length(naCoeff) > 0)
        {
            warning("covariates do not meet regularity conditions; ",
                    "the following covariates have been excluded: ",
                    paste(names(out@coefficients)[naCoeff], collapse=" "),
                    call.=FALSE)

            out@coefficients <- out@coefficients[-naCoeff]
            out@model.matrix <- moMat[, -naCoeff, drop=FALSE]
        }

        res <- try(solve(crossprod(out@model.matrix,
                                   out@model.matrix * out@variance)),
                   silent=TRUE)

        if (is(res, "try-error"))
            stop("covariates do not meet regularity conditions => ",
                 "possibly remove constant and/or linearly dependent columns",
                 call.=FALSE)
        else
            out@inv.matrix <- res

        if (adj == "force" || (adj == "automatic" &&  length(glmres$y) < 2000))
        {
            message("small sample correction applied")

            out@res.resampling.adj <- resampling(out, n.resampling.adj)

            if (adjExact)
            {
                VX <- out@model.matrix * out@variance
                P0 <- -VX %*% tcrossprod(out@inv.matrix, VX)
                diag(P0) <- diag(P0) + out@variance

                e <- eigen(P0, symmetric=TRUE)

                out@P0sqrt <- tcrossprod(sweep(e$vectors, MARGIN=2,
                                               STATS=sqrt(pmax(e$values, 0)),
                                               FUN="*"),
                                         e$vectors)
            }
        }

        if (n.resampling > 0)
            out@res.resampling <- resampling(out, n.resampling)
    }
    else if (type == "linear")
    {
        lmres <- lm(X, data=y)

        moMat <- model.matrix(lmres)

        out <- new("NullModel",
                   type=type,
                   residuals=lmres$residuals,
                   model.matrix=moMat,
                   coefficients=lmres$coefficients,
                   na.omit=omitted,
                   n.cases=0,
                   type.resampling=type.resampling,
                   prob=numeric(),
                   variance=summary(lmres)$sigma^2)

        naCoeff <- which(is.na(out@coefficients))

        if (length(naCoeff) > 0)
        {
            warning("covariates do not meet regularity conditions; ",
                    "the following covariates have been excluded: ",
                    paste(names(out@coefficients)[naCoeff], collapse=" "),
                    call.=FALSE)

            out@coefficients <- out@coefficients[-naCoeff]
            out@model.matrix <- moMat[, -naCoeff, drop=FALSE]
        }

        res <- try(solve(crossprod(moMat)), silent=TRUE)

        if (is(res, "try-error"))
            stop("covariates do not meet regularity conditions => ",
                 "possibly remove constant and/or linearly dependent columns",
                 call.=FALSE)
        else
            out@inv.matrix <- res

        if (n.resampling > 0)
            out@res.resampling <- resampling(out, n.resampling)
    }
    else ## type == "bernoulli"
    {
        y <- mf[, 1]
        names(y) <- rownames(mf)
        p <- mean(y)

        if (n.resampling > 0)
            warning("resampling not available for type 'exact' => ",
                    "ignoring 'n.resampling' argument", call.=FALSE)

        out <- new("NullModel",
                   type=type,
                   residuals=y,
                   na.omit=omitted,
                   n.cases=sum(y),
                   prob=p,
                   variance=(p * (1 - p)))
    }

    if (length(names(out@residuals)) == 0)
        warning("no sample names available => null model ",
                "will not be usable for association test with VCF file",
                call.=FALSE)

    out@call <- deparse(sys.call(-1))

    out
}

setMethod("nullModel", signature(X="formula", y="missing"), nullModel.formula)
setMethod("nullModel", signature(X="formula", y="data.frame"),
          nullModel.formula)
