summary.VariantInfo <- function(object)
{
    cat("Variant info:\n")
    cat("\tNumber of variants:", length(object), "\n\n")

    res <- list(meanMAF=NA, medianMAF=NA, minMAF=NA, maxMAF=NA, tableTypes=NULL)

    if (length(object) > 0)
    {
        MAF <- MAF(object)

        res$meanMAF   <- mean(MAF, na.rm=TRUE)
        res$medianMAF <- median(MAF, na.rm=TRUE)
        res$minMAF    <- min(MAF, na.rm=TRUE)
        res$maxMAF    <- max(MAF, na.rm=TRUE)

        cat("\tMean MAF:   ", res$meanMAF,   "\n")
        cat("\tMedian MAF: ", res$medianMAF, "\n")
        cat("\tMinimum MAF:", res$minMAF,    "\n")
        cat("\tMaximum MAF:", res$maxMAF,    "\n")

        if (length(mcols(object)$type) > 0)
        {
            res$tableTypes <- table(mcols(object)$type, useNA="no")
            n <- sum(res$tableTypes)

            if (n > 0)
            {
                w1 <- max(nchar(names(res$tableTypes)))
                w2 <- ceiling(log10(max(res$tableTypes)))

                labels <- format(paste0(names(res$tableTypes), ":"),
                                 width=(w1 + 2))
                absNo  <- format(res$tableTypes, justify="right", width=w2)
                perc   <- format(res$tableTypes * 100 / n, digits=3)

                cat("\n",
                    paste("\t", labels, absNo, " (", perc, "%)",
                          sep="",collapse="\n"), sep="", "\n")
            }
            else
                cat("\n\tonly NAs in metadata column 'type'\n")
        }
        else
            cat("\n\tno metadata column 'type' available\n")
    }

    invisible(res)
}

setMethod("summary", signature(object="VariantInfo"), summary.VariantInfo)
