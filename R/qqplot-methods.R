qqplot.vector <- function(x, y, xlab="", ylab="", xlim=NA, ylim=NA,
                          common.scale=TRUE, lwd=1, lcol="red", ...,
                          conf.level=NULL, conf.args=NULL)
{
    if (missing(y))
    {
        y <- sort(x)
        x <- ppoints(length(y))
    }
    else
    {
        if (length(x) != length(y))
            stop("'x' and 'y' contain different numbers of regions",
                 call.=FALSE)

        x <- sort(x)
        y <- sort(y)
    }

    if (common.scale)
    {
        xlim <- c(0, max(c(-log10(x[which(x > 0)]),
                           -log10(y[which(y > 0)]))))
        ylim <- xlim
    }
    else
    {
        if (identical(xlim, NA))
            xlim <- c(0, max(-log10(x[which(x > 0)])))

        if (identical(ylim, NA))
            ylim <- c(0, max(-log10(y[which(y > 0)])))
    }

    plot(NULL, NULL, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)

    if (lwd > 0)
        abline(0, 1, col=lcol, lwd=lwd)

    points(-log10(x), -log10(y), ...)

    invisible(list(x=x, y=y))
}


qqplot.AssocTestResultRanges <- function(x, y, xlab=deparse(substitute(x)),
                                         ylab=deparse(substitute(y)),
                                         common.scale=TRUE,
                                         preserveLabels=FALSE,
                                         lwd=1, lcol="red", ...,
                                         conf.level=NULL, conf.args=NULL)
{
    if (missing(y))
    {
        if (!preserveLabels)
        {
            if (is.null(xlab) || is.na(xlab))
            {
                xlab <- expression(paste("Expected ", ~~ -log[10](italic(p))))
                ylab <- expression(paste("Observed ", ~~ -log[10](italic(p))))
            }
            else if (is.character(xlab) || is.expression(xlab))
            {
                if (nchar(xlab) > 0)
                {
                    if (!preserveLabels)
                    {
                        ylab <- bquote(paste(-log[10](italic(p)), " in ",
                                             .(xlab)))
                        xlab=expression(paste("Expected ", -log[10](italic(p))))
                    }
                }
                else
                {
                    xlab <- ""
                    ylab <- ""
                }
            }
        }

        subsel <- which(mcols(x)$n > 0)

        qqplot.vector(mcols(x)$p.value[subsel], xlab=xlab, ylab=ylab,
                      common.scale=common.scale, lwd=lwd, lcol=lcol, ...)
    }
    else
    {
        if (length(x) != length(y))
            stop("'x' and 'y' contain different numbers of regions",
                 call.=FALSE)

        if (any(mcols(x)$n != mcols(y)$n))
            stop("numbers of variants in the regions of 'x' and 'y' ",
                 "do not match", call.=FALSE)

        if (!preserveLabels)
        {
            if (is.null(xlab) || is.na(xlab))
                xlab <- expression(paste("Observed ", -log[10](italic(p))))
            else if (is.character(xlab) || is.expression(xlab))
            {
                if (nchar(xlab) > 0)
                     if (!preserveLabels)
                        xlab <- bquote(paste(-log[10](italic(p)), " in ",
                                             .(xlab)))
                else
                    xlab <- ""
            }

            if (is.null(ylab) || is.na(ylab))
                ylab <- expression(paste("Observed ", -log[10](italic(p))))
            else if (is.character(ylab) || is.expression(ylab))
            {
                if (nchar(ylab) > 0)
                     if (!preserveLabels)
                        ylab <- bquote(paste(-log[10](italic(p)), " in ",
                                             .(ylab)))
                else
                    ylab <- ""
            }
        }

        subsel <- which(mcols(x)$n > 0)

        qqplot.vector(mcols(x)$p.value[subsel], mcols(y)$p.value[subsel],
                      xlab=xlab, ylab=ylab, common.scale=common.scale,
                      lwd=lwd, lcol=lcol, ...)
    }
}

setMethod("qqplot", signature(x="AssocTestResultRanges", y="missing"),
          qqplot.AssocTestResultRanges)

setMethod("qqplot",
          signature(x="AssocTestResultRanges", y="AssocTestResultRanges"),
          qqplot.AssocTestResultRanges)
