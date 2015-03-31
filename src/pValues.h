//
// File pValues.h defining prototypes for pValues.cpp
//

#ifndef _pValues_H_

#define _pValues_H_

RcppExport SEXP davies(SEXP lambdaR, SEXP xR, SEXP accR);
RcppExport SEXP liu(SEXP lambdaR, SEXP xR);
RcppExport SEXP liuMod(SEXP lambdaR, SEXP xR);

#endif
