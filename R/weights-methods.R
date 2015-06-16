weights.AssocTestResult <- function(object, Z, model)
{
    if (!is(Z, "GenotypeMatrix"))
        stop("'Z' must be a 'GenotypeMatrix' object", call.=FALSE)

    if (substr(object@kernel, 1, 6) != "linear")
        stop("variant contributions can only be computed for ",
             "kernels \"linear.podkat\" and \"linear.SKAT\"", call.=FALSE)

    if (!is(model, "NullModel"))
        stop("'model' must be a 'NullModel' object", call.=FALSE)

    if (!is.null(object@weights) && length(object@weights) != ncol(Z))
        stop("length of weight vector in 'object' does not match 'ncol(Z)'",
             call.=FALSE)

    if (length(object@samples) > 0 && length(object@samples) != length(model))
        stop("number of samples in 'object' does not fit to null model",
             call.=FALSE)

    ## check for correct order of samples
    if (length(rownames(Z)) > 0 && length(names(model@residuals)) > 0)
    {
        sel <- match(names(model@residuals), rownames(Z))

        if (any(is.na(sel)))
            stop("not all samples of null model contained in genotype matrix",
                 call.=FALSE)

        Z <- Z[sel, ]
    } ## if no names are present, assume that order is correct
    else if (nrow(Z) == (length(model) + length(model@na.omit)))
    {
        if (length(model@na.omit) > 0)
            Z <- Z[-model@na.omit, ]
    }
    else if (nrow(Z) != length(model))
        stop("dimension of genotype matrix 'Z' does not fit to null model",
             call.=FALSE)

    out <- variantInfo(Z)

    w <- computeWeights(Z, model@residuals, object@kernel,
                        object@weights, object@width)
    w2 <- w^2

    mcols(out) <- NULL
    mcols(out)$weight.raw <- w
    mcols(out)$weight.contribution <- w2 / sum(w2)

    class(out) <- "GRanges"

    out
}

setMethod("weights", signature(object="AssocTestResult"),
          weights.AssocTestResult)


weights.AssocTestResultRanges.GenotypeMatrix <- function(object, Z, model)
{
    n <- length(object)

    weights <- NULL

    if (is.numeric(object@weights))
    {
        if (length(object@weights) != ncol(Z))
            stop("length of weight vector in 'object' does not match ",
                 "'ncol(Z)'", call.=FALSE)

        weights <- object@weights
    }
    else if (is.function(object@weights))
    {
        if (length(MAF(Z)) > 0)
            weights <- object@weights(MAF(Z))
        else
            weights <- numeric()
    }

    ## check for correct order of samples
    if (length(rownames(Z)) > 0 && length(names(model@residuals)) > 0)
    {
        sel <- match(names(model@residuals), rownames(Z))

        if (any(is.na(sel)))
            stop("not all samples of null model contained in genotype matrix",
                 call.=FALSE)

        Z <- Z[sel, ]
    } ## if no names are present, assume that order is correct
    else if (nrow(Z) == (length(model) + length(model@na.omit)))
    {
        if (length(model@na.omit) > 0)
            Z <- Z[-model@na.omit, ]
    }
    else if (nrow(Z) != length(model))
        stop("dimension of genotype matrix 'Z' does not fit to null model",
             call.=FALSE)

    computeWeightsForRange <- function(i)
    {
        sel <- which(variantInfo(Z) %over% object[i])

        if (length(sel) > 0)
        {
            res <- variantInfo(Z)[sel]
            mcols(res) <- NULL

            w <- computeWeights(Z[, sel], residuals(model),
                                kernel=object@kernel, weights[sel],
                                object@width)
            w2 <- w^2

            mcols(res)$weight.raw <- w
            mcols(res)$weight.contribution <- w2 / sum(w2)

            class(res) <- "GRanges"
            res
        }
        else
            GRanges()
    }

    out <- GRangesList(lapply(1:n, computeWeightsForRange))
    seqlevels(out) <- seqlevels(object)
    seqlengths(out) <- seqlengths(object)
    genome(out) <- genome(object)

    names(out) <- paste0(seqnames(object), ":", start(object), "-", end(object))

    out
}

weights.AssocTestResultRanges.TabixFile <- function(object, Z, model, sex=NULL)
{
    ## check for data integrity and possible sub-setting
    sampleNames <- readSampleNamesFromVcfHeader(Z)

    if (length(sampleNames) == 0)
        stop("could not determine sample names from tabix file ",
             path(Z), call.=FALSE)

    if (length(object@samples) > 0 && length(object@samples) != length(model))
        stop("number of samples in 'object' does not fit to null model",
             call.=FALSE)
    else if (!all(object@samples %in% sampleNames))
        stop("some samples in 'object' do not occur in tabix file 'Z'",
             call.=FALSE)

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(names(sex)) == 0 && length(sex) != length(model))
            stop("length of 'sex' argument does not match number of ",
                 "samples in null model", call.=FALSE)
    }

    if (is.factor(sex))
    {
        lev <- levels(sex)

        if (length(lev) != 2 || any(lev != c("F", "M")))
            stop("invalid 'sex' argument", call.=FALSE)
    }

    ## try to determine sample names from null model
    sampleNamesNM <- names(model@residuals)

    subset <- logical()

    if (length(sampleNamesNM) > 0)
    {
        fileN <- length(sampleNames)

        subset <- (sampleNames %in% sampleNamesNM)

        if (length(which(subset)) != length(sampleNamesNM))
            stop("some samples of null model are missing in tabix file ",
                 path(Z), call.=FALSE)

        sampleNames <- sampleNames[subset]

        perm <- match(sampleNames, sampleNamesNM)

        model <- model[perm]

        if (!is.null(sex))
        {
            if (length(names(sex)) > 0)
            {
                perm <- match(sampleNames, names(sex))

                if (any(is.na(perm)))
                    stop("name mismatch between 'sex' and null model",
                         call.=FALSE)
            }

            sex <- sex[perm]
        }
    }
    else
        stop("no sample names in null model", call.=FALSE)

    Z <- readGenotypeMatrix(Z, reduce(object), which(subset),
                            noIndels=object@vcfParams$noIndels,
                            onlyPass=object@vcfParams$onlyPass,
                            na.limit=object@vcfParams$na.limit,
                            MAF.limit=object@vcfParams$MAF.limit,
                            na.action=object@vcfParams$na.action,
                            MAF.action=object@vcfParams$MAF.action,
                            sex=sex)

    weights(object, Z, model)
}


setMethod("weights", signature(object="AssocTestResultRanges"),
    function(object, Z, model, limit=20, sex=NULL)
    {
        if (length(object) == 0)
        {
            out <- GRangesList()
            mcols(out) <- list(weight.raw=numeric(),
                               weight.contribution=numeric())
            return(out)
        }

        if (substr(object@kernel, 1, 6) != "linear")
            stop("variant contributions can only be computed for ",
                 "kernels \"linear.podkat\" and \"linear.SKAT\"",
                 call.=FALSE)

        if (!is(model, "NullModel"))
            stop("'model' must be a 'NullModel' object", call.=FALSE)

        if (is.numeric(limit) && length(object) > limit)
            stop("more than ", limit, " regions in 'object'. ",
                 "Shorten 'object' or increase 'limit'.", call.=FALSE)

        if (is(Z, "GenotypeMatrix"))
            weights.AssocTestResultRanges.GenotypeMatrix(object, Z, model)
        else if (is(Z, "TabixFile"))
            weights.AssocTestResultRanges.TabixFile(object, Z, model, sex)
        else if (is(Z, "character"))
            weights.AssocTestResultRanges.TabixFile(object, TabixFile(Z),
                                                    model, sex)
        else
            stop("invalid 'Z' argument; must be 'GenotypeMatrix' ",
                 " object, 'TabixFile' object, or file name", call.=FALSE)
    })
