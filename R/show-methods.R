setMethod("show", signature(object="GenotypeMatrix"),
    function(object)
    {
        cat("Genotype matrix:\n")
        cat("\tNumber of samples:", nrow(object), "\n")
        cat("\tNumber of variants:", ncol(object), "\n\n")

        if (ncol(object) > 0)
        {
            MAF <- MAF(object)

            cat("\tMean MAF:   ", mean(MAF, na.rm=TRUE), "\n")
            cat("\tMedian MAF: ", median(MAF, na.rm=TRUE), "\n")
            cat("\tMinimum MAF:", min(MAF, na.rm=TRUE), "\n")
            cat("\tMaximum MAF:", max(MAF, na.rm=TRUE), "\n")
        }

        invisible(NULL)
    }
)

setMethod("show", signature(object="NullModel"),
    function(object)
    {
        if (object@type == "logistic")
            cat("Logistic model:\n")
        else if (object@type == "linear")
            cat("Linear model:\n")
        else
            cat("Simple Bernoulli model:\n")

        intercept <- as.numeric("(Intercept)" %in%
                                    colnames(object@model.matrix))

        ncov <- ncol(object@model.matrix) - intercept

        if (ncov == 0)
            if (intercept)
                cat("\tOnly intercept (no covariates)\n")
            else
                cat("\tRaw phenotypes (no covariates, no intercept)\n")
        else
            cat("\tNumber of covariates:", ncov,
                if (intercept) "(+ intercept)\n"
                else "(no intercept)\n")

        cat("\tNumber of samples:", length(object))

        if (length(object@na.omit) > 0)
            cat("\t(", length(object@na.omit), " samples omitted due to NAs)\n",
                sep="")
        else
            cat("\n")

        if (object@type == "logistic" || object@type == "bernoulli")
            cat("\tNumber of positives (cases):", object@n.cases, "\n")
        else
            cat("\tVariance of residuals:", object@variance, "\n")

        if (ncol(object@res.resampling) > 0)
            cat("\tResampling: ", ncol(object@res.resampling),
                " repeats (", object@type.resampling, ")\n", sep="")
        else
            cat("\tNo resampling\n")

        if (ncol(object@res.resampling.adj) > 0)
            cat("\tAdjustment of higher moments: ",
                ncol(object@res.resampling.adj),
                " repeats (", object@type.resampling, ")\n", sep="")
    }
)

setMethod("show", signature(object="AssocTestResult"),
    function(object)
    {
        cat("Association test results:\n")

        cat("\tNull model:", object@type, "\n")
        cat("\tNumber of samples:", object@dim[1], "\n")
        cat("\tNumber of variants:", object@dim[2], "\n")

        cat("\tKernel:", object@kernel, "\n")
        cat("\tTest statistic:", object@Q, "\n")

        cat("\tp-value:", object@p.value, "\n")

        if (object@type == "logistic" && length(object@correction) > 0)
        {
            if (all(object@correction))
                cat("\t(small sample correction +",
                    "correction for higher moments applied)\n")
            else if (object@correction["exact"])
                cat("\t(small sample correction applied)\n")
            else if (object@correction["resampling"])
                cat("\t(correction for higher moments applied)\n")
        }

        invisible(NULL)
    })

setMethod("show", signature(object="AssocTestResultRanges"),
          function(object)
              print(object, cutoff=NA, max.show=0))
