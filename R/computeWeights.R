## note: variants making use of sparseness have turned out to produce
##       additional overhead; that's why this function casts its input to
##       a regular dense matrix

computeWeights <- function(Z, res, kernel=c("linear.podkat", "linear.SKAT"),
                           weights=NULL, width=1000)
{
    kernel <- match.arg(kernel)

    pos <- start(Z@variantInfo)
    Z <- as.matrix(Z)

    kk <- unlist(strsplit(kernel, ".", fixed=TRUE))

    if (!is.null(weights))
        Z <- sweep(Z, MARGIN=2, STATS=weights, FUN="*")

    if (kernel == "linear.podkat")
        Z <- Z %*% .Call("posKernel", as.integer(pos), as.double(width),
                         PACKAGE="podkat")

    drop(t(res) %*% Z)
}
