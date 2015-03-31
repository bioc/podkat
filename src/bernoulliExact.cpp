#include <Rcpp.h>

using namespace Rcpp;

#include "bernoulliExact.h"


double traverseSummands(int n, double x, NumericMatrix &K, IntegerVector flag,
			double sum, double total, int pos, double p,
			double limit)
{
    int k;
    double travTotal = 0;
/*
    for (k = 0; k <= pos + 1; k++) Rprintf(" ");
    Rprintf("Position = %d / sum = %3.2f / total prob = %6.4g / Flags:\n",
	    pos, sum, total);

    for (k = 0; k <= pos + 1; k++) Rprintf(" ");
    for (k = 0; k < n; k++)
	Rprintf(" %d", flag[k]);

    Rprintf("\n");
*/
    for (k = pos + 1; k < n; k++)
    {
	int zeros = 1;
	flag[k] = 0;

	double nsum = sum;
	
	for (int j = 0; j < k; j++)
	{
	    if (flag[j])
		nsum -= (2 * K(k, j));
	    else
		zeros++;
	}

	nsum -= K(k, k);

	for (int j = k + 1; j < n; j++)
	    nsum -= (2 * K(k, j));

	if (nsum >= x)
	{
	    double local;

	    local = pow(p, n - zeros) * pow(1 - p, zeros);
	    local += traverseSummands(n, x, K, flag, nsum, total + local,
				      k, p, limit);

	    travTotal += local;

	    if (total + travTotal >= limit)
		return limit;
	}

	flag[k] = 1;
    }
/*    for (k = 0; k <= pos + 1; k++) Rprintf(" ");
      Rprintf("=> returning %6.4g\n", travTotal);*/
    return travTotal;
}


RcppExport SEXP computeExactBernoulliPvalue(SEXP xR, SEXP KR, SEXP pR,
					    SEXP limitR)
{
    NumericMatrix K(KR);
    double x = as<double>(xR), limit = as<double>(limitR), p = as<double>(pR);
    double sum = 0, total = 0;
    int n = K.nrow(), i, j;
    IntegerVector flag(n, 1);

    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    sum += K(i, j);

    total = traverseSummands(n, x, K, flag, sum, pow(p, n), -1, p, limit);

    return NumericVector::create(total);
}
