## note: variants making use of sparseness have turned out to produce
##       additional overhead; that's why this function casts its input to
##       a regular dense matrix

computeKernel <- function(Z, kernel=c("linear.podkat",
                                      "quadratic.podkat",
                                      "localsim.podkat",
                                      "linear.SKAT",
                                      "quadratic.SKAT",
                                      "localsim.SKAT"),
                          weights=NULL, pos=NULL, width=1000)
{
    kernel <- match.arg(kernel)

    if (is(Z, "GenotypeMatrix"))
    {
        if (is.null(weights))
        {
            if (length(unique(seqnames(Z@variantInfo))) > 1)
                stop("genotype matrix 'Z' contains variants from more than ",
                     "one chromosome", call.=FALSE)

            pos <- start(Z@variantInfo)
        }

        Z <- as.matrix(Z)
    }
    else if (is(Z, "Matrix"))
        Z <- as.matrix(Z)
    else if (!is(Z, "matrix") || !is.numeric(Z))
        stop("'Z' must be object of class 'matrix' or 'Matrix'", call.=FALSE)

    if (!is.null(weights))
    {
        if (is.numeric(weights))
        {
            if (length(weights) != ncol(Z))
                stop("length of 'weights' does not match 'ncol(Z)'",
                     call.=FALSE)
        }
        else
            stop("invalid 'weights' argument", call.=FALSE)
    }

    kk <- unlist(strsplit(kernel, ".", fixed=TRUE))

    if (kk[1] != "localsim" && !is.null(weights))
        Z <- sweep(Z, MARGIN=2, STATS=weights, FUN="*")

    if (kk[2] == "podkat")
    {
        if (is.null(pos))
            stop("positions not provided: cannot use kernel ", kernel,
                 call.=FALSE)
        else if (is.numeric(pos))
        {
            if (length(pos) != ncol(Z))
                stop("number of positions does not match number of variants",
                     call.=FALSE)
        }
        else
            stop("invalid 'pos' argument", call.=FALSE)

        if (!is.numeric(width) || length(width) != 1 || width <= 0)
            stop("invalid 'width' argument", call.=FALSE)

        Z <- Z %*% .Call("posKernel", as.integer(pos), as.double(width),
                         PACKAGE="podkat")
    }

    if (kk[1] == "linear")
        K <- tcrossprod(Z)
    else if (kk[1] == "quadratic")
        K <- (1 + tcrossprod(Z))^2
    else if (kk[1] == "localsim")
    {
        if (is.null(weights))
            K <- .Call("localSimKernel", Z, PACKAGE="podkat")
        else
            K <- .Call("localSimKernelWeighted", Z, as.double(weights^2),
                       PACKAGE="podkat")
    }

    attr(K, "kernel") < kernel

    K
}
