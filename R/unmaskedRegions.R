unmaskedRegions <- function(x, chrs=character(), pseudoautosomal=NULL,
                           ignoreGaps=250,
                           activeMasks=active(masks(x[[1]])))
{
    if (!is(x, "MaskedBSgenome"))
        stop("'x' must be 'MaskedBSgenome' object", call.=FALSE)

    if (!is.null(chrs) && !is.character(chrs))
        stop("invalid 'chrs' argument", call.=FALSE)

    if (!is.logical(activeMasks))
        stop("'activeMasks' must be a logical vector",
             call.=FALSE)

    if (length(chrs) > 0)
    {
        if (any(!(chrs %in% seqlevels(x))))
            warning("some chromosome names from 'chrs' were omitted that are ",
                    "not present in 'x'", call.=FALSE)

        chrs <- intersect(chrs, seqnames(x))
    }
    else
        chrs <- seqnames(x)

    pseudoauto <- GRanges()

    if (is.data.frame(pseudoautosomal))
    {
        colPresent <- c("chrom", "start.base", "end.base") %in%
                      colnames(pseudoautosomal)

        if (any(!colPresent))
            stop("not all required columns present in data frame ",
                 "'pseudoautosomal'", call.=FALSE)

        if ((!is.character(pseudoautosomal$chrom) &&
             !is.factor(pseudoautosomal$chrom)) ||
            !is.numeric(pseudoautosomal$start.base) ||
            !is.numeric(pseudoautosomal$end.base))
            stop("invalid column type in data frame 'pseudoautosomal'",
                 call.=FALSE)

        pseudoautosomal$chrom <- as.character(pseudoautosomal$chrom)

        if (!all(unique(pseudoautosomal$chrom) %in% seqlevels(x)))
            stop("invalid chromosome names in 'chrom' column of data frame ",
                 "'pseudoautosomal'", call.=FALSE)

        sel <- which(pseudoautosomal$chrom %in% chrs)

        if (length(sel) > 0)
        {
            pseudoauto <- GRanges(seqnames=pseudoautosomal$chrom,
                                  ranges=
                                  IRanges(start=pseudoautosomal$start.base,
                                          end=pseudoautosomal$end.base))

            names(pseudoauto) <- rownames(pseudoautosomal)
        }
        else
            pseudoauto <- GRanges()
    }
    else if (!is.null(pseudoautosomal) && !is.na(pseudoautosomal))
        stop("invalid 'pseudoautosomal' argument", call.=FALSE)

    unmasked <- GRangesList(lapply(chrs, function(chr)
                {
                    msk <- masks(x[[chr]])
                    active(msk) <- activeMasks
                    msk <- collapse(msk)

                    if (maskedwidth(msk) == 0)
                        GRanges(seqnames=chr,
                                ranges=IRanges(start=1, end=seqlengths(x)[chr]))
                    else
                    {
                        masked <- GRanges(seqnames=chr,
                                          ranges=IRanges(msk[[1, exact=TRUE]]))
                        ignore <- which((width(masked) < ignoreGaps) &
                                        (start(masked) != 1) &
                                        (end(masked) != seqlengths(x)[chr]))

                        if (length(ignore) > 0)
                        {
                            masked <- masked[-ignore]

                            if (length(masked) == 0)
                                return(GRanges(seqnames=chr,
                                               ranges=IRanges(start=1,
                                               end=seqlengths(x)[chr])))
                        }

                        gaps(masked)
                    }
                }))

    seqlengths(unmasked)[chrs] <- seqlengths(x)[chrs]

    names(unmasked) <- chrs

    if (length(pseudoauto) > 0)
    {
        for (i in 1:length(pseudoauto))
        {
            unmasked[[names(pseudoauto)[i]]] <-
                intersect(unmasked[[as.character(seqnames(pseudoauto)[i])]],
                          pseudoauto[i])

            unmasked[[as.character(seqnames(pseudoauto)[i])]] <-
                setdiff(unmasked[[as.character(seqnames(pseudoauto)[i])]],
                        pseudoauto[i])
        }

        names(unmasked) <- c(chrs, names(pseudoauto))
    }

    seqinfo(unmasked) <- seqinfo(x)[chrs]

    unmasked
}
