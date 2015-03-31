c.AssocTestResultRanges <- function(x, ..., recursive=FALSE)
{
    if (!identical(recursive, FALSE))
        stop("'recursive' argument not supported", call.=FALSE)

    if (missing(x))
        args <- unname(list(...))
    else
        args <- unname(list(x, ...))

    if (length(args) == 0)
        return(NULL)

    type <- unique(sapply(args, function(object) object@type))

    if (length(type) > 1)
        stop("results can only be merged if they are based on the same ",
             "type of null model", call.=FALSE)

    genome <- unique(sapply(args, function(object) genome(object)[1]))

    if (length(genome) > 1)
        stop("results can only be merged if they are based on the same genome",
             call.=FALSE)

    kernel <- unique(sapply(args, function(object) object@kernel))

    if (length(kernel) > 1)
        stop("results can only be merged if they have been computed using ",
             "the same kernel", call.=FALSE)

    width <- unique(sapply(args, function(object) object@width))

    if (length(width) > 1)
        stop("results can only be merged if they have been computed using ",
             "the same width parameter", call.=FALSE)

    samples1 <- args[[1]]@samples

    lapply(args[-1],
           function(object)
               if (length(object@samples) != length(samples1) ||
                   any(object@samples != samples1))
                   stop("results can only be merged if they have been ",
                        "computed for the same samples", call.=FALSE))

    wClass <- unique(sapply(args, function(object) class(object@weights)))

    if (wClass == "function")
        weights <- args[[1]]@weights
    else if (wClass == "numeric")
        weights <- do.call(c, lapply(args, function(object) object@weights))

    if (length(wClass) > 1)
        stop("results can only be merged if they have been computed using ",
             "the same weighting method", call.=FALSE)

    vcfParams <- lapply(args, function(object) object@vcfParams)

    if (length(unique(sapply(vcfParams, length))) > 1)
        stop("results can only be merged if they have been computed using ",
             "the same VCF file reading parameters", call.=FALSE)

    if (length(vcfParams[[1]]) > 0)
    {
        for (para in names(vcfParams[[1]]))
            if (length(unique(sapply(vcfParams,
                                     function(object) object[[para]]))) > 1)
                stop("results can only be merged if they have been computed ",
                     "using the same VCF file reading parameters", call.=FALSE)
    }

    sexN <- sapply(args, function(object) length(object@sex))

    if (length(which(sexN > 0)) > 0)
        sex <- args[[which(sexN > 0)[1]]]@sex
    else
        sex <- factor()

    args <- lapply(args, function(object)
                         {
                             mcols(object)$p.value.adj <- NULL
                             mcols(object)$p.value.resampled.adj <- NULL
                             as(object, "GRanges")
                         })

    x <- as(unlist(GRangesList(args)), "AssocTestResultRanges")

    x@type       <- type
    x@samples    <- samples1
    x@kernel     <- kernel
    x@weights    <- weights
    x@width      <- width
    x@vcfParams  <- vcfParams[[1]]
    x@sex        <- sex

    x@call       <- deparse(sys.call(-1))
    x@adj.method <- "none"

    x
}

setMethod("c", signature("AssocTestResultRanges"),
          c.AssocTestResultRanges)
