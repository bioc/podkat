readRegionsFromBedFile <- function(file, header=FALSE, sep="\t",
                                   col.names=c("chrom", "chromStart",
                                               "chromEnd", "names"),
                                   ignoreMcols=TRUE, seqInfo=NULL)
{
    df <- read.table(file=file, header=header, sep=sep,
                     stringsAsFactors=FALSE, col.names=col.names)

    if (length(df$chrom) == 0)
        stop("'chrom' column not present in BED file", call.=FALSE)

    if (length(df$chromStart) == 0)
        stop("'chromStart' column not present in BED file", call.=FALSE)

    if (length(df$chromEnd) == 0)
        stop("'chromEnd' column not present in BED file", call.=FALSE)

    gr <- GRanges(seqnames=df$chrom,
                  ranges=IRanges(start=df$chromStart, end=df$chromEnd))

    if ("strand" %in% col.names)
        strand(gr) <- df$strand

    if (length(df$names) > 0)
        names(gr) <- df$names

    if (!ignoreMcols)
    {
        mCol <- setdiff(colnames(df), c("chrom", "chromStart", "chromEnd",
                                        "names", "strand", "width"))

        if (length(mCol) > 0)
            mcols(gr) <- df[, mCol]
    }

    if (!is.null(seqInfo) && is(seqInfo, "Seqinfo"))
    {
        seqlevels(gr) <- seqlevels(seqInfo)
        seqlengths(gr) <- seqlengths(seqInfo)
        genome(gr) <- genome(seqInfo)
    }

    gr
}
