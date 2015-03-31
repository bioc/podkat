partitionRegions.GRanges <- function(x, chrs=character(), width=5000,
                                     overlap=0.5)
{
    overlap <- check.arg.numeric(overlap, 0, FALSE, 0.8, FALSE)
    width   <- check.arg.integer(width, 250, FALSE, NA, FALSE)

    if (width > 100000)
        warning("widths > 100,000 not recommended", call.=FALSE)

    if (!is.null(chrs) && !is.character(chrs))
        stop("invalid 'chrs' argument", call.=FALSE)

    if (length(chrs) > 0)
    {
        names(chrs) <- NULL

        if (any(!(chrs %in% seqlevels(x))))
            warning("some undefined chromosome names were omitted", call.=FALSE)

        x <- x[which(as.character(seqnames(x)) %in% chrs)]
    }

    if (length(x) == 0)
        segments <- GRanges(seqinfo=seqinfo(x))
    else
    {
        res <- .Call("partitionRegions", as.character(seqnames(x)),
                     as.integer(start(x)), as.integer(end(x)),
                     as.integer(width), as.double(overlap))

        segments <- GRanges(seqnames=res$seqnames,
                            ranges=IRanges(start=res$start, end=res$end),
                            seqinfo=seqinfo(x))
    }

    segments
}

setMethod("partitionRegions", signature=signature(x="GRanges"),
          partitionRegions.GRanges)


setMethod("partitionRegions", signature=signature(x="GRangesList"),
          function(x, chrs=character(), width=5000, overlap=0.5)
          {
              if (!is.null(chrs) && !is.character(chrs))
                  stop("invalid 'chrs' argument", call.=FALSE)

              if (length(chrs) > 0)
              {
                  names(chrs) <- NULL

                  if (any(!(chrs %in% seqlevels(x))))
                      warning("some undefined chromosome names were omitted",
                              call.=FALSE)
              }

              out <- GRangesList(lapply(x, partitionRegions.GRanges,
                                        chrs=chrs, width=width,
                                        overlap=overlap))

              mcols(out) <- mcols(x)

              out[which(sapply(out, length) > 0)]
          })


setMethod("partitionRegions", signature=signature(x="MaskedBSgenome"),
          function(x, chrs=character(), width=5000, overlap=0.5, ...)
              partitionRegions(unmaskedRegions(x, chrs, ...), width=width,
                               overlap=overlap))
