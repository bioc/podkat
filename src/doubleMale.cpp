#include <Rcpp.h>

using namespace Rcpp;

#include "doubleMale.h"


/*-------------------------------------------------------------------------*/
/* auxiliary function for multiplying male samples' genotypes with 2 in    */
/*     a sparse matrix in an efficient way                                 */
/*-------------------------------------------------------------------------*/

RcppExport SEXP doubleMale(SEXP iR, SEXP xR, SEXP sexR)
{
    DoubleVector x(xR);
    IntegerVector i(iR), sex(sexR);
    int I, N = x.length();
    DoubleVector newX(N);

    for (I = 0; I < N; I++)
    {
	if (sex[i[I]] > 1 && x[I] <= 1)
	    newX[I] = 2 * x[I];
	else
	    newX[I] = x[I];
    }

    return newX;
}
