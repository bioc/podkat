betaWeights <- function(x, shape1=1, shape2=25)
{
    if (missing(x))
        function(x) dbeta(x, shape1=shape1, shape2=shape2)
    else
        dbeta(x, shape1=shape1, shape2=shape2)
}

logisticWeights <- function(x, th=0.07, slope=150)
{
    if (missing(x))
        function(x) 1 / (1 + exp(slope * (x - th)))
    else
        1 / (1 + exp(slope * (x - th)))
}

invSdWeights <- function(x)
{
    if (missing(x))
        function(x) 1 / sqrt(x * (1 - x))
    else
        1 / sqrt(x * (1 - x))
}
