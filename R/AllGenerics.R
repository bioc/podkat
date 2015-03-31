setGeneric("assocTest", package="podkat",
           function(Z, model, ...) standardGeneric("assocTest"))

setGeneric("nullModel", package="podkat",
           function(X, y, ...) standardGeneric("nullModel"))

setGeneric("genotypeMatrix", package="podkat",
           function(Z, pos, seqnames, ...) standardGeneric("genotypeMatrix"))

setGeneric("variantInfo", package="podkat",
           function(x) standardGeneric("variantInfo"))

setGeneric("MAF", package="podkat",
           function(x) standardGeneric("MAF"))

setGeneric("partitionRegions", package="podkat",
           function(x, ...) standardGeneric( "partitionRegions"))

setGeneric("readGenotypeMatrix", function(file, regions, ...)
           standardGeneric("readGenotypeMatrix"))

setGeneric("readVariantInfo", function(file, regions, ...)
           standardGeneric("readVariantInfo"))

setGeneric("filterResult", function(object, ...)
           standardGeneric("filterResult"))
