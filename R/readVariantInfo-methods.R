.vcfScanVariantInfo <- function(file, seqnames, start, end, subset, noIndels,
                                onlyPass, na.limit, MAF.limit, na.action,
                                MAF.action, refAlt, omitZeroMAF, sex)
{
    if (!Rsamtools::isOpen(file))
    {
        open(file)
        on.exit(close(file))
    }

    result <- .Call("readVariantInfo", file$.extptr,
                    as.character(seqnames), as.integer(start),
                    as.integer(end), as.logical(subset),
                    as.logical(noIndels), as.logical(onlyPass),
                    as.double(na.limit), as.double(MAF.limit),
                    as.integer(switch(na.action,
                                      impute.major=1,
                                      omit=2,
                                      fail=4)),
                    as.integer(switch(MAF.action,
                                      invert=1,
                                      omit=2,
                                      ignore=3,
                                      fail=4)),
                    as.logical(refAlt), as.logical(omitZeroMAF),
                    as.integer(sex),
                    PACKAGE="podkat")

    if (is.null(result))
        return(variantInfo())
    else if (is.character(result))
        stop(result, call.=FALSE)
    else if (is.list(result))
    {
        vInfo <- GRanges(seqnames=result$seqname,
                         ranges=IRanges(start=result$pos, width=1))
        names(vInfo) <- result$names
        mcols(vInfo) <- data.frame(type=result$type, MAF=result$MAF)
        class(vInfo) <- "VariantInfo"

        if (refAlt)
        {
            mcols(vInfo)$ref <- result$ref
            mcols(vInfo)$alt <- result$alt
        }

        vInfo
    }
    else
        result
}

readVariantInfo.TabixFile <- function(file, regions, subset,
                                      noIndels=TRUE, onlyPass=TRUE,
                                      na.limit=1, MAF.limit=1,
                                      na.action=c("impute.major",
                                                  "omit", "fail"),
                                      MAF.action=c("ignore", "omit",
                                                   "invert", "fail"),
                                      omitZeroMAF=TRUE, refAlt=FALSE, sex=NULL)
{
    ## input checks
    na.action   <- match.arg(na.action)
    MAF.action  <- match.arg(MAF.action)
    noIndels    <- check.arg.logical(noIndels)
    onlyPass    <- check.arg.logical(onlyPass)
    omitZeroMAF <- check.arg.logical(omitZeroMAF)
    refAlt      <- check.arg.logical(refAlt)
    MAF.limit   <- check.arg.numeric(MAF.limit, 0, TRUE, 1, FALSE)
    na.limit    <- check.arg.numeric(na.limit, 0, TRUE, 1, FALSE)

    sampleNames <- readSampleNamesFromVcfHeader(file)

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (missing(subset) && length(sex) != length(sampleNames))
            stop("length of 'sex' argument does not match number of ",
                 "samples in null model", call.=FALSE)
    }

    if (is.factor(sex))
    {
        lev <- levels(sex)

        if (length(lev) != 2 || any(lev != c("F", "M")))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(names(sex)) > 0)
        {
            sel <- match(sampleNames, names(sex))

            if (any(is.na(sel)))
                stop("sample name mismatch between VCF file and ",
                     "'sex' argument", call.=FALSE)

            sex <- sex[sel]
        }
    }

    if (missing(subset) || is.null(subset))
        subset <- logical()
    else if (is.character(subset))
    {
        subsetT <- (sampleNames %in% subset)

        if (length(which(subsetT)) != length(subset))
            stop("some sample names in 'subset' are missing in VCF file ",
                 path(file), call.=FALSE)

        subset <- subsetT

        if (!is.null(sex))
            sex <- sex[subset]
    }
    else if (is.numeric(subset) || is.integer(subset))
    {
        subsetT <- (1:length(sampleNames) %in% subset)

        if (length(which(subsetT)) != length(subset))
            stop("some sample indices in 'subset' are out of range",
                 call.=FALSE)

        subset <- subsetT

        if (!is.null(sex))
            sex <- sex[subset]
    }
    else if (!is.null(subset))
        stop("invalid 'subset' argument", call.=FALSE)

    .vcfScanVariantInfo(file, seqnames=as.character(seqnames(regions)),
                        start=start(regions), end=end(regions),
                        subset=subset, noIndels=as.logical(noIndels),
                        onlyPass=as.logical(onlyPass),
                        MAF.limit=MAF.limit, na.limit=na.limit,
                        na.action=na.action, MAF.action=MAF.action,
                        refAlt=as.logical(refAlt),
                        omitZeroMAF=as.logical(omitZeroMAF), sex=sex)
}


setMethod("readVariantInfo", signature(file="TabixFile", regions="GRanges"),
    readVariantInfo.TabixFile)

setMethod("readVariantInfo", signature(file="TabixFile", regions="missing"),
    function(file, regions, ...) readVariantInfo(file, GRanges(), ...))

setMethod("readVariantInfo", signature(file="character", regions="missing"),
    function(file, regions, ...) readVariantInfo(TabixFile(file), GRanges(),
                                                 ...))
