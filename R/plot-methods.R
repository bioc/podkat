plot.AssocTestResultRanges.character <-
    function(x, y, cutoff=0.05,
             which=c("p.value", "p.value.adj", "p.value.resampled",
                     "p.value.resampled.adj"),
             showEmpty=FALSE, as.dots=FALSE, pch=19,
             col=c("darkgrey", "grey"), scol="red", lcol="red",
             xlab=NULL, ylab=NULL, ylim=NULL, lwd=1, cex=1,
             cexXaxs=1, cexYaxs=1, srt=0, adj=c(0.5, 1), ...)
{
    which <- match.arg(which)

    if (is.null(mcols(x)[[which]]))
        stop("column '", which, "' missing", call.=FALSE)

    if (!showEmpty)
        x <- x[which(mcols(x)$n > 0)]

    if (missing(y) || length(y) == 0)
        chrs <- seqlevels(x)
    else if (is.null(y) || is.na(y))
        chrs <- levels(factor(as.character(seqnames(x))))
    else
    {
        chrs <- y

        if (!any(chrs %in% seqlevels(x)))
            stop("'chrs' argument contains invalid sequence names", call.=FALSE)

        x <- x[which(as.character(seqnames(x)) %in% chrs)]
    }

    if (!is.null(ylim) && !is.na(ylim))
    {
        if (!is.numeric(ylim) || length(ylim) != 2 || ylim[1] != 0 ||
            ylim[2] <= 0)
            stop("ylim must be numeric vector of type c(0, x) with x > 0")
    }

    len <- seqlengths(x)[chrs]

    if (any(is.na(len)))
    {
        warning("some sequence lengths are not available in 'x'", call.=FALSE)

        correctSeqLength <- function(chr)
        {
            if (is.na(len[chr]))
                max(end(x)[which(seqnames(x) == chr)])
            else
                len[chr]
        }

        len <- sapply(chrs, correctSeqLength, USE.NAMES=FALSE)

        names(len) <- chrs
    }

    if (length(chrs) == 1)
    {
        if (is.null(xlab) || is.na(xlab))
            xlab <- paste("Chromosome", chrs, "of", genome(x)[1])

        rng <- GRanges(seqnames=chrs, ranges=IRanges(start=1, end=len))

        return(plot.AssocTestResultRanges.GRanges(x=x, y=rng, cutoff=cutoff,
                                                  which=which,
                                                  showEmpty=showEmpty,
                                                  as.dots=as.dots,
                                                  pch=pch, col=col[1],
                                                  scol=scol, lcol=lcol,
                                                  xlab=xlab, ylab=ylab,
                                                  ylim=ylim, lwd=lwd, cex=cex,
                                                  cexXaxs=cexXaxs,
                                                  cexYaxs=cexYaxs,
                                                  ...))
    }

    xlim <- c(1, sum(as.numeric(len)))

    pValues <- mcols(x)[[which]]

    if (is.null(ylim) || is.na(ylim))
        ylim <- c(0, max(-log10(pValues[which(pValues > 0)])))

    if (is.null(xlab) || is.na(xlab))
        xlab <- paste("Chromosomes of", genome(x)[1])

    if (is.null(ylab) || is.na(ylab))
        ylab <- expression(-log[10](italic(p)))

    endX <- cumsum(as.numeric(len))
    startX <- c(0, endX[-length(endX)]) + 1
    names(startX) <- chrs
    names(endX) <- chrs

    cols <- rep(col, length=length(chrs))
    names(cols) <- chrs

    plot(NULL, NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE,
         ...)

    tickMarks <- c(startX[1], (startX[-1] + endX[-length(endX)]) / 2,
                   endX[length(endX)])

    box()
    axis(1, at=tickMarks, labels=FALSE)
    axis(3, at=tickMarks, labels=FALSE)
    axis(2, cex.axis=cexYaxs)

    yC <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.03
    text(x=((startX + endX) / 2), y=yC, labels=chrs,
         srt=srt, adj=adj, xpd=TRUE, cex=cexXaxs)

    if (is.numeric(lwd) && lwd > 0 && is.numeric(cutoff) && cutoff > 0)
    {
        axis(4, at=-log10(cutoff), labels=FALSE, col.ticks=lcol,
             lwd.ticks=lwd, cex.axis=cexYaxs)
        mtext(format.pval(cutoff), 4, 1, at=-log10(cutoff), col=lcol,
              cex=cexYaxs)
        abline(-log10(cutoff), 0, lwd=lwd, col=lcol)
    }

    seqNames <- as.character(seqnames(x))

    plotChr <- function(chr)
    {
        sel <- which(seqNames == chr)

        colS <- ifelse(pValues[sel] <= cutoff, scol, cols[chr])

        if (as.dots)
            points(startX[chr] + (start(x)[sel] + end(x)[sel]) / 2,
                   -log10(pValues[sel]), pch=pch, col=colS, cex=cex, ...)
        else
            segments(x0=(startX[chr] + start(x)[sel]),
                     y0=-log10(pValues[sel]), x1=(startX[chr] + end(x)[sel]),
                     col=colS, lwd=lwd, ...)
    }

    ret <- sapply(chrs, plotChr)

    return(invisible(ylim))
}

setMethod("plot", signature(x="AssocTestResultRanges", y="missing"),
          plot.AssocTestResultRanges.character)

setMethod("plot", signature(x="AssocTestResultRanges", y="character"),
          plot.AssocTestResultRanges.character)


plot.AssocTestResultRanges.GRanges <-
    function(x, y, cutoff=0.05,
             which=c("p.value", "p.value.adj", "p.value.resampled",
                     "p.value.resampled.adj"),
             showEmpty=FALSE, as.dots=FALSE, pch=19,
             col="darkgrey", scol="red", lcol="red",
             xlab=NULL, ylab=NULL, ylim=NULL, lwd=1, cex=1,
             cexXaxs=1, cexYaxs=1, ...)
{
    which <- match.arg(which)

    if (is.null(mcols(x)[[which]]))
        stop("column '", which, "' missing", call.=FALSE)

    if (length(y) != 1)
        stop("only single genomic region can be plotted", call.=FALSE)

    if (!is.null(ylim) && !is.na(ylim))
    {
        if (!is.numeric(ylim) || length(ylim) != 2 || ylim[1] != 0 ||
            ylim[2] <= 0)
            stop("'ylim' must be numeric vector of type 'c(0, x)' with x > 0",
                 call.=FALSE)
    }

    x <- x[x %over% y]

    if (!showEmpty)
        x <- x[mcols(x)$n > 0]

    xlim <- c(start(y), end(y))

    pValues <- mcols(x)[[which]]

    if (is.null(ylim) || is.na(ylim))
        ylim <- c(0, max(-log10(pValues[which(pValues > 0)])))

    if (is.null(xlab) || is.na(xlab))
    {
        if (length(names(y)) > 0)
            xlab <- names(y)
        else
            xlab <- paste("Region [", start(y), "..", end(y),
                          "] of chromosome ", seqnames(y), " in ", genome(x)[1],
                          sep="")
    }

    if (is.null(ylab) || is.na(ylab))
        ylab <- expression(-log[10](italic(p)))

    plot(NULL, NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE,
         ...)

    tickMarks <- axTicks(1)

    box()
    axis(1, at=tickMarks, labels=format(as.integer(tickMarks), big.mark=",",
                                        trim=TRUE), cex.axis=cexXaxs)
    axis(3, at=xlim, labels=format(as.integer(xlim), big.mark=",",
                                   trim=TRUE), cex.axis=cexXaxs)
    axis(2, cex.axis=cexYaxs)

    if (is.numeric(lwd) && lwd > 0 && is.numeric(cutoff) && cutoff > 0)
    {
        axis(4, at=-log10(cutoff), labels=FALSE, col.ticks=lcol,
             lwd.ticks=lwd)
        mtext(format.pval(cutoff), 4, 1, at=-log10(cutoff), col=lcol,
              cex=cexYaxs)
        abline(-log10(cutoff), 0, lwd=lwd, col=lcol)
    }

    colS <- ifelse(pValues <= cutoff, scol, col)

    if (as.dots)
        points((start(x) + end(x)) / 2, -log10(pValues), pch=pch, col=colS,
               cex=cex, ...)
    else
        segments(x0=start(x), y0=-log10(pValues), x1=end(x), col=colS,
                 lwd=(10 * cex), ...)

    return(invisible(ylim))
}

setMethod("plot", signature(x="AssocTestResultRanges", y="GRanges"),
          plot.AssocTestResultRanges.GRanges)


plot.GenotypeMatrix.missing <- function(x, y, col="black",
                                        labRow=NULL, labCol=NULL,
                                        cexXaxs=(0.2 + 1 / log10(ncol(x))),
                                        cexYaxs=(0.2 + 1 / log10(nrow(x))),
                                        srt=90, adj=c(1, 0.5))
{
    sv <- variantInfo(x)

    x <- as.matrix(x)

    cv <- col2rgb(col[1])

    plot(NULL, NULL, xlim=c(1, ncol(x)), ylim=c(1, nrow(x)),
         axes=FALSE, xlab="", ylab="")

    image(x=1:ncol(x), y=1:nrow(x), z=t(x), zlim=c(0, 2),
          col=c("white", rgb(cv[1], cv[2], cv[3], 127, maxColorValue=255),
                rgb(cv[1], cv[2], cv[3], maxColorValue=255)),
          axes=FALSE, add=TRUE)

    if (is.null(labRow) && length(rownames(x)) > 0)
        labRow <- rownames(x)

    if (length(labRow) > 0)
        axis(2, at=1:nrow(x), labels=labRow, las=2, line=-0.5, tick=FALSE,
             cex.axis=cexYaxs)

    if (is.null(labCol) && length(names(sv)) > 0)
        labCol <- names(sv)

    if (length(labCol) > 0)
    {
        yC <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.03
        text(x=1:length(labCol), y=yC, labels=labCol,
             srt=srt, adj=adj, xpd=TRUE, cex=cexXaxs)
    }
}

setMethod("plot", signature(x="GenotypeMatrix", y="missing"),
          plot.GenotypeMatrix.missing)


plot.GenotypeMatrix.numeric <- function(x, y, col="black", ccol="red", lwd=2,
                                        labRow=NULL, labCol=NULL,
                                        cexXaxs=(0.2 + 1 / log10(ncol(x))),
                                        cexYaxs=(0.2 + 1 / log10(nrow(x))),
                                        srt=90, adj=c(1, 0.5))
{
    sv <- variantInfo(x)

    if (length(rownames(x)) > 0 && length(names(y)) > 0)
    {
        sel <- match(names(y), rownames(x))

        if (any(is.na(sel)))
            stop("'y' contains samples that are missing from 'x'", call.=FALSE)

        x <- as.matrix(x)[sel, ]
    }
    else if (length(y) != nrow(x))
        stop("'nrow(x)' does not match length of 'y'", call.=FALSE)
    else
        x <- as.matrix(x)

    perm <- order(y)
    x <- x[perm, ]
    y <- y[perm]

    plot(NULL, NULL, xlim=c(1, ncol(x)), ylim=c(1, nrow(x)),
         axes=FALSE, xlab="", ylab="")

    cv <- col2rgb(col[1])

    image(x=1:ncol(x), y=1:nrow(x), z=t(x), zlim=c(0, 2),
          col=c("white", rgb(cv[1], cv[2], cv[3], 127, maxColorValue=255),
                rgb(cv[1], cv[2], cv[3], maxColorValue=255)),
          axes=FALSE, add=TRUE)

    if (is.null(labRow) && length(rownames(x)) > 0)
        labRow <- rownames(x)

    if (length(labRow) > 0)
        axis(2, at=1:nrow(x), labels=labRow, las=2, line=-0.5, tick=FALSE,
             cex.axis=cexYaxs)

    if (is.null(labCol) && length(names(sv)) > 0)
        labCol <- names(sv)

    if (length(labCol) > 0)
    {
        yC <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.03
        text(x=1:length(labCol), y=yC, labels=labCol,
             srt=srt, adj=adj, xpd=TRUE, cex=cexXaxs)
    }

    par(new=TRUE)

    plot(y, 1:nrow(x), col=ccol, lwd=lwd, type="l", axes=FALSE,
         xlab="", ylab="")
    axis(3)
}

setMethod("plot", signature(x="GenotypeMatrix", y="numeric"),
          plot.GenotypeMatrix.numeric)


plot.GenotypeMatrix.factor <- function(x, y, col=rainbow(length(levels(y))),
                                       labRow=NULL, labCol=NULL,
                                       cexXaxs=(0.2 + 1 / log10(ncol(x))),
                                       cexYaxs=(0.2 + 1 / log10(nrow(x))),
                                       srt=90, adj=c(1, 0.5))
{
    sv <- variantInfo(x)

    if (length(rownames(x)) > 0 && length(names(y)) > 0)
    {
        sel <- match(names(y), rownames(x))

        if (any(is.na(sel)))
            stop("'y' contains samples that are missing from 'x'", call.=FALSE)

        x <- as.matrix(x)[sel, ]
    }
    else if (length(y) != nrow(x))
        stop("'nrow(x)' does not match length of 'y'", call.=FALSE)
    else
        x <- as.matrix(x)

    perm <- order(y)
    x <- x[perm ,]
    y <- y[perm]

    if (length(col) < length(levels(y)))
        col <- rep(col, length.out=length(levels(y)))

    plot(NULL, NULL, xlim=c(1, ncol(x)), ylim=c(1, nrow(x)),
         axes=FALSE, xlab="", ylab="")

    pos <- 0

    for (i in 1:length(levels(y)))
    {
        cv <- col2rgb(col[i])

        sel <- which(y == levels(y)[i])

        image(x=1:ncol(x), y=(pos + (1:length(sel))), z=t(x[sel, ]),
              zlim=c(0, 2),
              col=c("white", rgb(cv[1], cv[2], cv[3], 127, maxColorValue=255),
                    rgb(cv[1], cv[2], cv[3], maxColorValue=255)),
              axes=FALSE, add=TRUE)

        pos <- pos + length(sel)
    }

    if (is.null(labRow) && length(rownames(x)) > 0)
        labRow <- rownames(x)

    if (length(labRow) > 0)
        axis(2, at=1:nrow(x), labels=labRow, las=2, line=-0.5, tick=FALSE,
             cex.axis=cexYaxs)

    if (is.null(labCol) && length(names(sv)) > 0)
        labCol <- names(sv)

    if (length(labCol) > 0)
    {
        yC <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.03
        text(x=1:length(labCol), y=yC, labels=labCol,
             srt=srt, adj=adj, xpd=TRUE, cex=cexXaxs)
    }
}

setMethod("plot", signature(x="GenotypeMatrix", y="factor"),
          plot.GenotypeMatrix.factor)


plot.GRanges.character <- function(x, y, alongGenome=FALSE,
                                   type=c("r", "s", "S", "l", "p",
                                          "b", "c", "h", "n"),
                                   xlab=NULL, ylab=NULL, col="red", lwd=2,
                                   cexXaxs=(0.2 + 1 / log10(length(x))),
                                   cexYaxs=1, frame.plot=TRUE,
                                   srt=90, adj=c(1, 0.5), ...)
{
    type <- match.arg(type)

    if (length(y) != 1 || !(y %in% colnames(mcols(x))))
        stop("'y' must be a single string referring to a metadata ",
             "column of 'x'")

    colData <- mcols(x)[[y]]

    if (!is.numeric(colData))
        stop("column '", y, "' is not numeric")

    if (is.null(ylab) || is.na(ylab))
        ylab <- y

    if (alongGenome)
    {
        if (length(unique(seqnames(x))) != 1)
            stop("plotting along genome only possible if all regions are from ",
                 "the same sequence")

        if (is.null(xlab) || is.na(xlab))
            xlab <- seqnames(x)[1]

        if (type == "r")
        {
            plot(NULL, NULL, xlab=xlab, ylab=ylab,
                 xlim=c(min(start(x)), max(end(x))),
                 ylim=range(colData, na.rm=TRUE), frame.plot=frame.plot,
                 axes=FALSE, ...)
            segments(x0=start(x), y0=colData, x1=end(x), col=col, lwd=lwd)

            axis(1, cex.axis=cexXaxs)
            axis(2, cex.axis=cexYaxs)
        }
        else
        {
            pos <- (start(x) + end(x)) / 2
            sel <- order(pos)

            plot(pos[sel], colData[sel], type=type, xlab=xlab, ylab=ylab,
                 col=col, lwd=lwd, frame.plot=frame.plot, axes=FALSE, ...)

            axis(1, cex.axis=cexXaxs)
            axis(2, cex.axis=cexYaxs)
        }
    }
    else
    {
        if (is.null(xlab) || is.na(xlab))
            xlab <- ""

        if (type == "r")
        {
            plot(NULL, NULL, axes=FALSE, frame.plot=frame.plot,
                 xlab=xlab, ylab=ylab, xlim=c(0.5, length(x) + 0.5),
                 ylim=range(colData, na.rm=TRUE), ...)
            segments(x0=1:length(x) - 0.5, y0=colData, x1=1:length(x) + 0.5,
                     col=col, lwd=lwd)
        }
        else
        {
            plot(1:length(x), colData, type=type, xlab=xlab, ylab=ylab,
                 axes=FALSE, frame.plot=frame.plot, col=col, lwd=lwd, ...)
        }

        axis(2, cex.axis=cexYaxs)

        if (length(names(x)) > 0)
        {
            yC <- par("usr")[3] - (par("usr")[4] - par("usr")[3]) * 0.03
            text(x=1:length(x), y=yC, labels=names(x),
                 srt=srt, adj=adj, xpd=TRUE, cex=cexXaxs)
        }
        else
            axis(1, cex.axis=cexXaxs)
    }
}

setMethod("plot", signature(x="GRanges", y="character"),
          plot.GRanges.character)
