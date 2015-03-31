readRegionsFromBedFile <- function(file, header=FALSE, sep="\t",
                                   col.names=c("chrom", "chromStart",
                                               "chromEnd", "names", "width",
                                               "strand"),
                                   seqInfo=NULL)
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

    if (length(df$names) > 0)
        names(gr) <- df$names

    if (!is.null(seqInfo) && is(seqInfo, "Seqinfo"))
    {
        seqlevels(gr) <- seqlevels(seqInfo)
        seqlengths(gr) <- seqlengths(seqInfo)
        genome(gr) <- genome(seqInfo)
    }

    gr
}
