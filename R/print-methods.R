print.AssocTestResultRanges <- function(x, cutoff=0.05,
                                        sortBy=c("p.value", "p.value.adj",
                                                 "p.value.resampled",
                                                 "p.value.resampled.adj",
                                                 "genome", "none"),
                                        max.show=10)
{
    sortBy <- match.arg(sortBy)

    if (substr(sortBy, 1, 7) == "p.value" && is.null(mcols(x)[[sortBy]]))
        stop("column '", sortBy, "' missing", call.=FALSE)

    cat("Overview of association test:\n")

    cat("\tNull model:", x@type, "\n")
    cat("\tNumber of samples:", length(x@samples), "\n")
    cat("\tNumber of regions:", length(x), "\n")

    empty <- which(mcols(x)$n == 0)

    cat("\tNumber of regions without variants:", length(empty), "\n")
    cat("\tAverage number of variants in regions:",
        format(mean(mcols(x)$n), digits=1, nsmall=1),
        "\n")
    cat("\tGenome:", genome(x)[1], "\n")
    cat("\tKernel:", x@kernel, "\n")
    cat("\tp-value adjustment:", x@adj.method, "\n")

    if (is.numeric(cutoff) && !any(is.na(cutoff)) && length(cutoff) > 0)
    {
        cat("\nOverview of significance of results:\n")

        cutoff <- sort(cutoff, decreasing=FALSE)

        ret <- sapply(cutoff,
                      function(th)
                      cat("\tNumber of tests with p < ", th, ": ",
                          length(which(mcols(x)$p.value < th)),
                          "\n", sep=""))

        if (length(mcols(x)$p.value.adj) > 0)
            ret <- sapply(cutoff,
                          function(th)
                          cat("\tNumber of tests with adj. p < ", th, ": ",
                              length(which(mcols(x)$p.value.adj < th)),
                              "\n", sep=""))
    }

    if (is.numeric(max.show) && max.show > 0)
    {
        if (sortBy == "genome")
        {
            rnk <- order(as(x, "GRanges"))
            x <- x[rnk]

            cat("\nResults for the first", min(length(x), max.show),
                "regions along the genome:\n")
        }
        else if (sortBy != "none")
        {
            x <- x[which(mcols(x)$n > 0)]
            if (is.numeric(cutoff) && length(cutoff) >= 1)
                x <- x[which(mcols(x)[[sortBy]] <= max(cutoff))]
            x <- x[order(mcols(x)[[sortBy]])]

            cat("\nResults for the", min(length(x), max.show),
                "most significant regions:\n")
        }
        else
            cat("\nResults for the first", min(length(x), max.show),
                "regions:\n")

        df <- as(x, "data.frame")

        selCols <- c("seqnames", "start", "end", "width",
                     colnames(mcols(x)))

        show(head(df[, selCols], max.show))
    }

    invisible(x)
}

setMethod("print", signature(x="AssocTestResultRanges"),
          print.AssocTestResultRanges)
