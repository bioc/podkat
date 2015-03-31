filterResult.AssocTestResultRanges <-
    function(object, cutoff=0.05,
             filterBy=c("p.value", "p.value.adj", "p.value.resampled",
                        "p.value.resampled.adj"))
{
    filterBy <- match.arg(filterBy)

    if (is.null(mcols(object)[[filterBy]]))
        stop("column '", filterBy, "' missing", call.=FALSE)

    object[which(mcols(object)[[filterBy]] <= cutoff)]
}

setMethod("filterResult", signature(object="AssocTestResultRanges"),
          filterResult.AssocTestResultRanges)


filterResult.GRanges <- function(object, cutoff=0.1)
{
    if (length(mcols(object)$weight.contribution) == 0)
        stop("'object' has no metadata column 'weight.contribution'",
                call.=FALSE)
    else
       object[which(mcols(object)$weight.contribution >= cutoff)]
}

setMethod("filterResult", signature(object="GRanges"),
          filterResult.GRanges)

setMethod("filterResult", signature(object="GRangesList"),
          function(object, cutoff=0.1)
              GRangesList(lapply(object, filterResult.GRanges, cutoff=cutoff)))
