genotypeMatrix.dgCMatrix.GRanges <- function(Z, pos, ploidy=2,
                                             na.limit=1, MAF.limit=1,
                                             na.action=c("impute.major",
                                                         "omit", "fail"),
                                             MAF.action=c("invert", "omit",
                                                          "ignore", "fail"),
                                             sex=NULL)
{
    na.action   <- match.arg(na.action)
    MAF.action  <- match.arg(MAF.action)
    MAF.limit   <- check.arg.numeric(MAF.limit, 0, TRUE, 1, FALSE)
    na.limit    <- check.arg.numeric(na.limit, 0, TRUE, 1, FALSE)
    ploidy      <- check.arg.integer(ploidy, 1, FALSE, NA, FALSE)

    if (length(pos) != ncol(Z))
        stop("'length(pos)' must equal 'ncol(Z)'", call.=FALSE)

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(sex) != nrow(Z))
            stop("length of 'sex' argument does not match number of ",
                 "rows of 'Z'", call.=FALSE)
    }

    if (is.factor(sex))
    {
        lev <- levels(sex)

        if (length(lev) != 2 || any(lev != c("F", "M")))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(names(sex)) > 0 && length(rownames(Z)) > 0)
        {
            sel <- match(sampleNames, names(sex))

            if (any(is.na(sel)))
                stop("sample name mismatch between VCF object and ",
                     "'sex' argument", call.=FALSE)

            sex <- sex[sel]
        }
    }

    ## multiply males' genotypes with 2 if user wants
    if (!is.null(sex) && ploidy == 2)
        Z@x <- .Call("doubleMale", as.integer(Z@i), as.double(Z@x),
                     as.integer(sex), PACKAGE="podkat")

    result <- .Call("checkAndFixGenotype",
                    as.integer(Z@i), as.integer(Z@p), as.integer(Z@Dim),
                    as.double(Z@x), as.integer(ploidy),
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
                    PACKAGE="podkat")

    if (is.character(result))
        stop(result, call.=FALSE)

    out <- new("dgCMatrix",
               i=result$i,
               p=result$p,
               x=result$x,
               Dim=result$Dim)

    sel <- which(!result$omitted)

    if (length(rownames(Z)) > 0)
        rownames(out) <- rownames(Z)

    if (length(names(pos)) > 0)
    {
        if (length(colnames(Z)) > 0 && !all(names(pos) == colnames(Z)))
            warning("column names of 'Z' are not consistent with names of ",
                    "'pos'; giving priority to the names of 'pos'", call.=FALSE)

        colnames(out) <- names(pos)[sel]
    }
    else if (length(colnames(Z)) > 0)
        colnames(out) <- colnames(Z)[sel]

    vInfo <- pos[sel]
    names(vInfo) <- NULL
    mcols(vInfo)$MAF <- result$MAF
    class(vInfo) <- "VariantInfo"

    out <- as(out, "GenotypeMatrix")
    out@variantInfo <- vInfo

    out
}

genotypeMatrix.Cmatrix.GRanges <- function(Z, pos,
                                           na.string=NULL,
                                           na.limit=1, MAF.limit=1,
                                           na.action=c("impute.major", "omit",
                                                       "fail"),
                                           MAF.action=c("invert", "omit",
                                                        "ignore", "fail"),
                                           sex=NULL)
{
    na.action   <- match.arg(na.action)
    MAF.action  <- match.arg(MAF.action)
    MAF.limit   <- check.arg.numeric(MAF.limit, 0, TRUE, 1, FALSE)
    na.limit    <- check.arg.numeric(na.limit, 0, TRUE, 1, FALSE)

    if (length(pos) != ncol(Z))
        stop("'length(pos)' must equal 'ncol(Z)'", call.=FALSE)

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(sex) != nrow(Z))
            stop("length of 'sex' argument does not match number of ",
                 "rows of 'Z'", call.=FALSE)
    }

    if (is.factor(sex))
    {
        lev <- levels(sex)

        if (length(lev) != 2 || any(lev != c("F", "M")))
            stop("invalid 'sex' argument", call.=FALSE)
    }

    if (!is.null(na.string))
    {
        if (is.character(na.string))
            Z[grep(".", Z, fixed=TRUE)] <- na.string
        else
            stop("invalid 'na.string' argument", call.=FALSE)
    }

    result <- .Call("checkAndFixGenotypeChar",
                    Z, as.double(na.limit), as.double(MAF.limit),
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

    if (is.character(result))
        stop(result, call.=FALSE)

    out <- new("dgCMatrix",
               i=result$i,
               p=result$p,
               x=result$x,
               Dim=result$Dim)

    sel <- which(!result$omitted)

    if (length(rownames(Z)) > 0)
        rownames(out) <- rownames(Z)

    if (length(names(pos)) > 0)
    {
        if (length(colnames(Z)) > 0 && !all(names(pos) == colnames(Z)))
            warning("column names of 'Z' are not consistent with names of ",
                    "'pos'; giving priority to the names of 'pos'", call.=FALSE)

        colnames(out) <- names(pos)[sel]
    }
    else if (length(colnames(Z)) > 0)
        colnames(out) <- colnames(Z)[sel]

    vInfo <- pos[sel]
    names(vInfo) <- NULL
    mcols(vInfo)$MAF <- result$MAF
    class(vInfo) <- "VariantInfo"

    out <- as(out, "GenotypeMatrix")
    out@variantInfo <- vInfo

    out
}

setMethod("genotypeMatrix",
          signature(Z="ANY", pos="GRanges", seqnames="missing"),
          function(Z, pos, seqnames, ploidy=2, na.string=NULL,
                   na.limit=1, MAF.limit=1,
                   na.action=c("impute.major", "omit", "fail"),
                   MAF.action=c("invert", "omit", "ignore", "fail"), sex=NULL)
          {
              if (is(Z, "dgCMatrix"))
                  genotypeMatrix.dgCMatrix.GRanges(Z=Z, pos=pos, ploidy=ploidy,
                                                   na.limit=na.limit,
                                                   MAF.limit=MAF.limit,
                                                   na.action=na.action,
                                                   MAF.action=MAF.action,
                                                   sex=sex)
              else if (is.matrix(Z) && is.numeric(Z))
                  genotypeMatrix.dgCMatrix.GRanges(Z=as(Z, "dgCMatrix"),
                                                   pos=pos, ploidy=ploidy,
                                                   na.limit=na.limit,
                                                   MAF.limit=MAF.limit,
                                                   na.action=na.action,
                                                   MAF.action=MAF.action,
                                                   sex=sex)
              else if (is.matrix(Z) && is.character(Z))
                  genotypeMatrix.Cmatrix.GRanges(Z=Z, pos=pos,
                                                 na.string=na.string,
                                                 na.limit=na.limit,
                                                 MAF.limit=MAF.limit,
                                                 na.action=na.action,
                                                 MAF.action=MAF.action,
                                                 sex=sex)
              else
                  stop("'Z' must be 'dgCMatrix' object or a matrix of mode ",
                       "'numeric' or 'character'", call.=FALSE)
          })

setMethod("genotypeMatrix",
          signature(Z="ANY", pos="numeric", seqnames="character"),
          function(Z, pos, seqnames, ...)
          {
              pos <- as.integer(pos)

              if (length(seqnames) != length(pos) && length(seqnames) != 1)
                  stop("'length(seqnames)' must be 1 or equal 'length(pos)'",
                       call.=FALSE)

              gr <- GRanges(seqnames=seqnames,
                            ranges=IRanges(start=pos, width=1))

              genotypeMatrix(Z, pos=gr, ...)
          })


setMethod("genotypeMatrix",
          signature(Z="ANY", pos="character", seqnames="missing"),
          function(Z, pos, seqnames, ...)
          {
              raw <- strsplit(pos, ":", fixed=TRUE)

              if (any(sapply(raw, length) != 2))
                  stop("invalid 'pos' argument, follow 'chr:pos' syntax",
                       call.=FALSE)

              raw <- unlist(raw)
              seqnames <- raw[seq(from=1, to=length(raw), by=2)]
              pos <- as.integer(raw[seq(from=2, to=length(raw), by=2)])

              gr <- GRanges(seqnames=seqnames,
                            ranges=IRanges(start=pos, width=1))

              genotypeMatrix(Z, pos=gr, ...)
          })


setMethod("genotypeMatrix",
          signature(Z="eSet", pos="character", seqnames="character"),
          function(Z, pos, seqnames, ...)
          {
              if (length(pos) != 1 || length(seqnames) != 1)
                  stop("both 'pos' and 'seqnames' must be strings identifying ",
                       "the feature data columns containing position data",
                       call.=FALSE)
              else if (length(featureData(Z)[[seqnames]]) == 0)
                  stop("no column '", seqnames, "' in 'featureData(Z)'",
                       call.=FALSE)
              else if (length(featureData(Z)[[pos]]) == 0)
                  stop("no column '", pos, "' in 'featureData(Z)'",
                       call.=FALSE)
              else if (is.null(assayData(Z)$call))
                  stop("no component 'call' in 'assayData(Z)'",
                       call.=FALSE)
              else
                  genotypeMatrix(as(t(assayData(Z)$call) - 1, "dgCMatrix"),
                                 as.numeric(featureData(Z)[[pos]]),
                                 as.character(featureData(Z)[[seqnames]]),
                                 ...)
          })


setMethod("genotypeMatrix",
          signature(Z="eSet", pos="numeric", seqnames="character"),
          function(Z, pos, seqnames, ...)
          {
              if (is.null(assayData(Z)$call))
                  stop("no component 'call' in 'assayData(Z)'",
                       call.=FALSE)
              else if (!is.matrix(assayData(Z)$call))
                  stop("'assayData(Z)$cal' is not a matrix", call.=FALSE)
              else
                  genotypeMatrix(as(t(assayData(Z)$call) - 1, "dgCMatrix"),
                                 pos, seqnames, ...)
          })


setMethod("genotypeMatrix",
          signature(Z="eSet", pos="character", seqnames="missing"),
          function(Z, pos, seqnames, ...)
          {
              if (is.null(assayData(Z)$call))
                  stop("no component 'call' in 'assayData(Z)'",
                       call.=FALSE)
              else if (!is.matrix(assayData(Z)$call))
                  stop("'assayData(Z)$cal' is not a matrix", call.=FALSE)
              else
                  genotypeMatrix(as(t(assayData(Z)$call) - 1, "dgCMatrix"),
                                 pos, ...)
          })


genotypeMatrix.ANY <- function(Z, pos, seqnames, subset,
                               noIndels=TRUE, onlyPass=TRUE, sex=NULL, ...)
{
    if (!is(Z, "VCF"))
        stop("if 'pos' and 'seqnames' are missing, 'Z' must be an object of ",
             "class 'VCF' (package 'VariantAnnotation')", call.=FALSE)
    else if (!requireNamespace("VariantAnnotation", quietly=TRUE))
        stop("could not load package 'VariantAnnotation'")

    if (length(VariantAnnotation::geno(Z)) == 0 ||
        length(VariantAnnotation::geno(Z)$GT) == 0)
        stop("genotype not included in VCF object", call.=FALSE)

    noIndels <- check.arg.logical(noIndels)
    onlyPass <- check.arg.logical(onlyPass)

    vSubset <- NULL

    if (onlyPass)
    {
        if (length(VariantAnnotation::fixed(Z)) == 0 ||
            length(VariantAnnotation::fixed(Z)$FILTER) == 0)
            warning("'FILTER' column not present in VCF object => ",
                    "ignoring 'onlyPass' argument", call.=FALSE)
        else
            vSubset <- which(VariantAnnotation::fixed(Z)$FILTER == "PASS")
    }

    if (noIndels)
    {
        if (length(VariantAnnotation::fixed(Z)) == 0 ||
            length(VariantAnnotation::fixed(Z)$REF) == 0)
            warning("'REF' column not present in VCF object => ",
                    "ignoring 'noIndels' argument", call.=FALSE)
        else
        {
            alt <- VariantAnnotation::fixed(Z)$ALT

            if (length(alt) == 0)
                warning("'ALT' column not present in VCF object => ",
                        "ignoring 'noIndels' argument", call.=FALSE)
            else
            {
                altLengths <- .Call("cumMax",
                                    as.integer(alt@unlistData@ranges@width),
                                    as.integer(alt@partitioning@end),
                                    PACKAGE="podkat")

                refLengths <- nchar(VariantAnnotation::fixed(Z)$REF)

                if (is.null(vSubset))
                    vSubset <- which(refLengths == 1 & altLengths == 1)
                else
                    vSubset <- intersect(vSubset,
                                         which(refLengths == 1 &
                                               altLengths == 1))
            }
        }
    }

    sampleNames <- samples(VariantAnnotation::header(Z))

    if (!is.null(sex))
    {
        if (is.character(sex))
            sex <- factor(sex)
        else if (!is.factor(sex))
            stop("invalid 'sex' argument", call.=FALSE)

        if (length(sex) != length(sampleNames))
            stop("length of 'sex' argument does not match number of ",
                 "samples in VCF object", call.=FALSE)
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
                stop("sample name mismatch between VCF object and ",
                     "'sex' argument", call.=FALSE)

            sex <- sex[sel]
        }
    }

    if (!missing(subset))
    {
        if (is.character(subset))
        {
            subset <- match(subset, sampleNames)

            if (any(is.na(sel)))
                stop("some sample names in subset argument are ",
                     "missing in VCF object", call.=FALSE)
        }
        else if (!is.numeric(subset) && !is.integer(subset))
            stop("invalid 'subset' argument", call.=FALSE)

        subset <- as.integer(subset)

        if (min(subset) < 1 || max(subset) > length(sampleNames))
            stop("some sample indices in 'subset' argument are out of range",
                 call.=FALSE)

        if (!is.null(sex))
            sex <- sex[subset]
    }
    else
        subset <- NULL

    if (is.null(vSubset))
    {
        if (is.null(subset))
            mat <- t(VariantAnnotation::geno(Z)$GT)
        else
            mat <- t(VariantAnnotation::geno(Z)$GT[, subset, drop=FALSE])

        gr <- GRanges(seqnames=rowRanges(Z)@seqnames,
                      ranges=rowRanges(Z)@ranges)
    }
    else
    {
        if (is.null(subset))
            mat <- t(VariantAnnotation::geno(Z)$GT[vSubset, , drop=FALSE])
        else
            mat <- t(VariantAnnotation::geno(Z)$GT[vSubset, subset, drop=FALSE])

        gr <- GRanges(seqnames=rowRanges(Z)@seqnames[vSubset],
                      ranges=rowRanges(Z)@ranges[vSubset])
    }

    seqlevels(gr) <- seqlevels(rowRanges(Z))
    seqlengths(gr) <- seqlengths(rowRanges(Z))
    genome(gr) <- genome(rowRanges(Z))

    genotypeMatrix.Cmatrix.GRanges(mat, gr, ...)
}

setMethod("genotypeMatrix",
          signature(Z="ANY", pos="missing", seqnames="missing"),
          genotypeMatrix.ANY)
