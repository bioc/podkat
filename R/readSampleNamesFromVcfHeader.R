readSampleNamesFromVcfHeader <- function(file, ...)
{
    if (is(file, "TabixFile"))
        file <- path(file)
    else if (!is(file, "character") || length(file) != 1)
        stop("'file' must either be a single path or a 'TabixFile' object",
             call.=FALSE)

    suppressWarnings(out <- scanBcfHeader(file, ...)[[1]]$Sample)

    out
}
