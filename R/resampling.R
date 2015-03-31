resampling <- function(nm, n)
{
    if (nm@type.resampling == "permutation")
        return(sapply(1:n, function(i) sample(nm@residuals)))
    else if (nm@type.resampling == "bootstrap")
    {
        l <- length(nm)

        if (nm@type == "linear")
            return(matrix(rnorm(l * n, mean=0, sd=sqrt(nm@variance)),
                          ncol=n))
        else
        {
            createRes <- function(i)
            {
                y1 <- rbinom(l, 1, nm@prob)
                y2 <- rbinom(l, 1, nm@prob)

                ident <- which(y1 & y2)
                sdiff <- which(y1 != y2)
                bzero <- which(!(y1 | y2))

                if (nm@n.cases <= length(ident))
                    samples <- sample(ident, nm@n.cases)
                else if (nm@n.cases <= length(ident) + length(sdiff))
                    samples <- c(ident,
                                 sample(sdiff, nm@n.cases - length(ident)))
                else
                    samples <- c(ident, sdiff,
                                 sample(bzero, nm@n.cases - length(ident)
                                                          - length(sdiff)))
                y <- -nm@prob
                y[samples] <- y[samples] + 1

                y
            }

            return(sapply(1:n, createRes))
        }
    }
    else
        stop("invalid resampling method", call.=FALSE)
}
