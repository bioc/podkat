setClass("NullModel",
    representation = representation
    (
        type               = "character",
        residuals          = "numeric",
        model.matrix       = "matrix",
        inv.matrix         = "matrix",
        P0sqrt             = "matrix",
        coefficients       = "numeric",
        na.omit            = "integer",
        n.cases            = "numeric",
        variance           = "numeric",
        prob               = "numeric",
        type.resampling    = "character",
        res.resampling     = "matrix",
        res.resampling.adj = "matrix",
        call               = "character"
    ),
    prototype = prototype
    (
        type               = "",
        residuals          = numeric(),
        model.matrix       = matrix(nrow=0, ncol=0),
        inv.matrix         = matrix(nrow=0, ncol=0),
        P0sqrt             = matrix(nrow=0, ncol=0),
        coefficients       = numeric(),
        na.omit            = integer(),
        n.cases            = 0,
        variance           = numeric(),
        prob               = numeric(),
        type.resampling    = "",
        res.resampling     = matrix(nrow=0, ncol=0),
        res.resampling.adj = matrix(nrow=0, ncol=0),
        call               = ""
    )
)

setClass("VariantInfo", contains="GRanges")

setClass("GenotypeMatrix",
    representation = representation
    (
        variantInfo = "VariantInfo"
    ),
    contains="dgCMatrix"
)

setClass("AssocTestResult",
    representation = representation
    (
        type               = "character",
        samples            = "character",
        kernel             = "character",
        dim                = "numeric",
        weights            = "numeric",
        width              = "numeric",
        method             = "ANY",
        correction         = "logical",
        Q                  = "numeric",
        p.value            = "numeric",
        Q.resampling       = "numeric",
        p.value.resampling = "numeric",
        p.value.resampled  = "numeric",
        call               = "character"
    ),
    prototype = prototype
    (
        type               = "",
        samples            = "",
        kernel             = "",
        dim                = c(0, 0),
        weights            = numeric(),
        width              = 0,
        method             = NULL,
        correction         = c(exact=FALSE, resampling=FALSE),
        Q                  = 0,
        p.value            = 0,
        Q.resampling       = 0,
        p.value.resampling = 0,
        p.value.resampled  = 0,
        call               = ""
     )
)

setClass("AssocTestResultRanges",
    representation = representation
    (
        type       = "character",
        samples    = "character",
        kernel     = "character",
        weights    = "ANY",
        width      = "numeric",
        adj.method = "character",
        vcfParams  = "list",
        sex        = "factor",
        call       = "character"
   ),
    prototype = prototype
    (
        type       = "",
        samples    = character(),
        kernel     = "",
        weights    = NULL,
        width      = 0,
        adj.method = "none",
        vcfParams  = list(),
        sex        = factor(),
        call       = ""
    ),
    contains="GRanges"
)
