p.adjust.AssocTestResultRanges <- function(p, method=p.adjust.methods,
                                           n=length(p))
{
    method <- match.arg(method)

    if (method == "none")
        return(p)

    p@adj.method <- method

    sel <- which(mcols(p)$n > 0)
    pval <- mcols(p)$p.value[sel]

    mcols(p)$p.value.adj <- 1

    mcols(p)$p.value.adj[sel] <- p.adjust(pval, method=method)

    if (length(mcols(p)$p.value.resampled) > 0)
    {
        sel <- which(mcols(p)$n > 0 &
                     !is.na(mcols(p)$p.value.resampled))
        pval <- mcols(p)$p.value.resampled[sel]

        mcols(p)$p.value.resampled.adj <- 1

        mcols(p)$p.value.resampled.adj[sel] <-
            p.adjust(pval, method=method, n=n)
    }

    p
}

setMethod("p.adjust", signature("AssocTestResultRanges"),
          p.adjust.AssocTestResultRanges)
