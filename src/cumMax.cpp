#include <Rcpp.h>

using namespace Rcpp;

#include "cumMax.h"

/*-------------------------------------------------------------------------*/
/* auxiliary function for determining the lengths of ALT alleles from the  */
/*    ALT column of VCF object (class 'DNAStringSetList')                  */
/*-------------------------------------------------------------------------*/

RcppExport SEXP cumMax(SEXP xR, SEXP pR)
{
    IntegerVector x(xR), p(pR);
    int i, index, maxi, n = p.length();
    IntegerVector res(n);

    for (i = 0, index = 0; i < n; i++)
    {
	maxi = -1;

	for (; index < p[i]; index++)
	    if (x[index] > maxi)
		maxi = x[index];

	res[i] = maxi;
    }

    return res;
}
