.vcfScan <- function(file, seqnames, start, end, subset, noIndels, onlyPass,
                     na.limit, MAF.limit, na.action, MAF.action, sex)
{
    if (!Rsamtools::isOpen(file))
    {
        open(file)
        on.exit(close(file))
    }

    result <- .Call("readGenotypeMatrix", file$.extptr,
                    as.character(seqnames), as.integer(start), as.integer(end),
                    as.logical(subset), as.logical(noIndels),
                    as.logical(onlyPass),
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
                    as.integer(sex),
                    PACKAGE="podkat")

    if (is.null(result))
        NULL
    else if (is.character(result))
        stop(result, call.=FALSE)
    else if (is.list(result))
    {
        GT <- new("dgCMatrix",
                  i=result$i,
                  p=result$p,
                  x=result$x,
                  Dim=result$Dim);

        list(seqnames=result$seqname, pos=result$pos, MAF=result$MAF,
             names=result$names, GT=GT)
    }
    else
        result
}

setMethod("readGenotypeMatrix", signature(file="TabixFile", regions="GRanges"),
    function(file, regions, subset, noIndels=TRUE, onlyPass=TRUE,
             na.limit=1, MAF.limit=1,
             na.action=c("impute.major", "omit", "fail"),
             MAF.action=c("invert", "omit", "ignore", "fail"), sex=NULL)
    {
        na.action <- match.arg(na.action)
        MAF.action <- match.arg(MAF.action)

        if (!is.numeric(MAF.limit) || length(MAF.limit) != 1 || MAF.limit > 1 ||
            MAF.limit <= 0)
            stop("'MAF.limit' must be numeric value above 0 and not ",
                 "larger than 1", call.=FALSE)

        if (!is.numeric(na.limit) || length(na.limit) != 1 || na.limit > 1 ||
            na.limit <= 0)
            stop("'na.limit' must be numeric value above 0 and not larger ",
                 "than 1", call.=FALSE)

        sampleNames <- readSampleNamesFromVcfHeader(file)

        if (!is.null(sex))
        {
            if (is.character(sex))
                sex <- factor(sex)
            else if (!is.factor(sex))
                stop("invalid 'sex' argument", call.=FALSE)

            if (missing(subset) && length(sex) != length(sampleNames))
                stop("length of 'sex' argument does not match number of ",
                     "samples in VCF file", call.=FALSE)
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
                stop("some sample names in 'subset' are ",
                     "missing in VCF file ", path(file),
                     call.=FALSE)

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

        result <- .vcfScan(file, seqnames=as.character(seqnames(regions)),
                           start=start(regions), end=end(regions),
                           subset=subset, noIndels=as.logical(noIndels),
                           onlyPass=as.logical(onlyPass),
                           MAF.limit=MAF.limit, na.limit=na.limit,
                           na.action=na.action, MAF.action=MAF.action, sex=sex)

        if (is.null(result))
        {
            out <-  as(as(matrix(nrow=max(length(subset), length(sampleNames)),
                                 ncol=0),
                          "dgCMatrix"),
                       "GenotypeMatrix")
        }
        else if (is.character(result))
            stop(result, call.=FALSE)
        else
        {
            out <- as(result$GT, "GenotypeMatrix")

            vInfo <- GRanges(seqnames=result$seqnames,
                             ranges=IRanges(start=result$pos, width=1))
            mcols(vInfo)$MAF <- result$MAF
            class(vInfo) <- "VariantInfo"

            out@variantInfo <- vInfo

            if (length(result$names) > 0)
                out@Dimnames[[2]] <- result$names
        }

        if (length(subset) > 0)
            out@Dimnames[[1]] <- sampleNames[subset]
        else
            out@Dimnames[[1]] <- sampleNames

        out
     })

setMethod("readGenotypeMatrix", signature(file="TabixFile", regions="missing"),
    function(file, regions, ...) readGenotypeMatrix(file, GRanges(), ...))

setMethod("readGenotypeMatrix", signature(file="character", regions="missing"),
    function(file, regions, ...) readGenotypeMatrix(TabixFile(file), ...))

setMethod("readGenotypeMatrix", signature(file="character", regions="GRanges"),
    function(file, regions, ...)
        readGenotypeMatrix(TabixFile(file), regions=regions, ...))
