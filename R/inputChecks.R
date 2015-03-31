check.arg.logical <- function(arg)
{
    argname <- deparse(substitute(arg))

    if (!is.logical(arg) || length(arg) != 1)
        stop("'", argname, "' must be single logical", call.=FALSE)

    arg
}

check.arg.numeric <- function(arg, lower=NA, strictLower=FALSE,
                              upper=NA, strictUpper=FALSE)
{
    argname <- deparse(substitute(arg))

    if (!is.numeric(arg) || length(arg) != 1)
        stop("'", argname, "' must be single numeric value", call.=FALSE)

    if (!is.na(lower))
    {
        if (strictLower)
        {
            if (arg <= lower)
                stop("'", argname, "' must be larger than ", lower, call.=FALSE)
        }
        else if (arg < lower)
            stop("'", argname, "' must be at least ", lower, call.=FALSE)
    }

    if (!is.na(upper))
    {
        if (strictUpper)
        {
            if (arg >= upper)
                stop("'", argname, "' must be smaller than ", upper,
                     call.=FALSE)
        }
        else if (arg > upper)
            stop("'", argname, "' must be at most ", upper, call.=FALSE)
    }

    arg
}

check.arg.integer <- function(arg, lower=NA, strictLower=FALSE,
                              upper=NA, strictUpper=FALSE)
{
    argname <- deparse(substitute(arg))

    if (!is.numeric(arg) || length(arg) != 1 || round(arg) != arg)
        stop("'", argname, "' must be single whole number", call.=FALSE)

    if (!is.na(lower))
    {
        if (strictLower)
        {
            if (arg <= lower)
                stop("'", argname, "' must be larger than ", lower, call.=FALSE)
        }
        else if (arg < lower)
            stop("'", argname, "' must be at least ", lower, call.=FALSE)
    }

    if (!is.na(upper))
    {
        if (strictUpper)
        {
            if (arg >= upper)
                stop("'", argname, "' must be smaller than ", upper,
                     call.=FALSE)
        }
        else if (arg > upper)
            stop("'", argname, "' must be at most ", upper, call.=FALSE)
    }

    arg
}
