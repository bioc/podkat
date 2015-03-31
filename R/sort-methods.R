setMethod("sort", signature(x="AssocTestResultRanges"),
    function(x, decreasing=FALSE,
             sortBy=c("p.value", "p.value.adj", "p.value.resampled",
                      "p.value.resampled.adj", "genome"))
    {
        sortBy <- match.arg(sortBy)

        if (sortBy == "genome")
            rnk <- order(as(x), "GRanges", decreasing=decreasing)
        else
        {
            if (is.null(mcols(x)[[sortBy]]))
                stop("column '", sortBy, "' missing", call.=FALSE)

            rnk <- order(mcols(x)[[sortBy]], decreasing=decreasing)
        }

        x[rnk]
    })
