#include <limits.h>
#include <float.h>
#include <Rcpp.h>

using namespace Rcpp;

#include "kernels.h"


RcppExport SEXP localSimKernel(SEXP ZR)
{
    NumericMatrix Z(ZR);
    int i, j, k, n = Z.nrow(), p = Z.ncol();
    NumericMatrix K(n, n);

    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    double val = 0;

	    for (k = 0; k < p; k++)
	    {
		double diff = 2 - fabs(Z(i, k) - Z(j, k));

		if (diff > 0)
		    val += diff;
	    }

	    K(i, j) = val / (2 * p);
	    K(j, i) = K(i, j);
	}

	K(i, i) = 1;
    }

    return(K);
}

RcppExport SEXP localSimKernelWeighted(SEXP ZR, SEXP weightsR)
{
    NumericMatrix Z(ZR);
    NumericVector weights(weightsR);
    int i, j, k, n = Z.nrow(), p = Z.ncol();
    NumericMatrix K(n, n);
    double wTotal = 0;

    for (k = 0; k < p; k++)
	wTotal += weights[k];

    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    double val = 0;

	    for (k = 0; k < p; k++)
	    {
		double diff = 2 - fabs(Z(i, k) - Z(j, k));

		if (diff > 0)
		    val += (weights[k] * diff);
	    }

	    K(i, j) = val / (2 * wTotal);
	    K(j, i) = K(i, j);
	}

	K(i, i) = 1;
    }

    return(K);
}

RcppExport SEXP posKernel(SEXP posR, SEXP widthR)
{
    NumericVector pos(posR);
    int n = pos.length(), i, j;
    double width = as<double>(widthR);
    NumericMatrix K(n, n);

    for (i = 0; i < n; i++)
    {
	for (j = (i + 1); j < n; j++)
	{
	    double dist = fabs(pos[i] - pos[j]);

	    if (dist < width)
	    {
		K(i, j) = 1 - dist / width;
		K(j, i) = K(i, j);
	    }
	    else
		break;
	}

	K(i, i) = 1;
    }

    return(K);
}
