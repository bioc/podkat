#include <limits.h>
#include <float.h>
#include <Rcpp.h>

#define SQRT2 1.414213562373095145475

using namespace Rcpp;

#include "pValues.h"
#include "qfc.h"


RcppExport SEXP davies(SEXP lambdaR, SEXP xR, SEXP accR)
{
    NumericVector lambda(lambdaR);
    NumericVector x(xR);
    double acc = as<double>(accR);
    int i, n = lambda.length(), m = x.length(), ifault, limit = 10000;
    NumericVector nc(n);
    IntegerVector deg(n, 1);
    NumericVector res(m);
    double c, result, trace[7], sigma = 0;
    
    for (i = 0; i < m; i++)
    {
	c = x[i];
	qfc(lambda.begin(), nc.begin(), deg.begin(), &n, &sigma, &c, &limit,
	    &acc, trace, &ifault, &result);

	if (ifault)
	    res[i] = -1;
	else
	    res[i] = 1. - result;
    }

    return(res);
}


RcppExport SEXP liu(SEXP lambdaR, SEXP xR)
{
    NumericVector lambda(lambdaR);
    NumericVector x(xR);
    int i, n = lambda.length(), m = x.length();
    NumericVector xArg(m);
    double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    double lsq, s1, s2, muQ, muX, sigmaQ, sigmaX, a, delta, l;

    for (i = 0; i < n; i++)
    {
	lsq = lambda[i] * lambda[i];

	c1 += lambda[i];
	c2 += lsq;
	c3 += lsq * lambda[i];
	c4 += lsq * lsq;
    }

    s1 = c3 / pow(c2, 1.5);
    s2 = c4 / (c2 * c2);

    muQ = c1;
    sigmaQ = sqrt(2 * c2);

    if ((s1 * s1) > s2)
    {
	a = 1 / (s1 - sqrt(s1 * s1 - s2));
	delta = a * a * (s1 * a - 1);
	l = a * a - 2 * delta;
    }
    else
    {
	a = 1 / s1;
	delta = 0;
	l = (c2 * c2 * c2) / (c3 * c3);
    }

    muX = l + delta;
    sigmaX = SQRT2 * a;

    for (i = 0; i < m; i++)
	xArg[i] = (x[i] - muQ) * sigmaX / sigmaQ + muX;

    return(wrap(pchisq(xArg, l, false, false)));
}


RcppExport SEXP liuMod(SEXP lambdaR, SEXP xR)
{
    NumericVector lambda(lambdaR);
    NumericVector x(xR);
    int i, n = lambda.length(), m = x.length();
    NumericVector xArg(m);
    double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    double lsq, s1, s2, muQ, muX, sigmaQ, sigmaX, a, delta, l;

    for (i = 0; i < n; i++)
    {
	lsq = lambda[i] * lambda[i];

	c1 += lambda[i];
	c2 += lsq;
	c3 += lsq * lambda[i];
	c4 += lsq * lsq;
    }

    s1 = c3 / pow(c2, 1.5);
    s2 = c4 / (c2 * c2);

    muQ = c1;
    sigmaQ = sqrt(2 * c2);

    if ((s1 * s1) > s2)
    {
	a = 1 / (s1 - sqrt(s1 * s1 - s2));
	delta = a * a * (s1 * a - 1);
	l = a * a - 2 * delta;
    }
    else
    {
	delta = 0;
	l = 1 / s2;
	a = sqrt(l);
    }

    muX = l + delta;
    sigmaX = SQRT2 * a;

    for (i = 0; i < m; i++)
	xArg[i] = (x[i] - muQ) * sigmaX / sigmaQ + muX;

    return(wrap(pchisq(xArg, l, false, false)));
}
