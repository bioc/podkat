//
// File kernels.h defining prototypes for kernels.cpp
//

#ifndef _kernels_H_

#define _kernels_H_

RcppExport SEXP localSimKernel(SEXP ZR);
RcppExport SEXP localSimKernelWeighted(SEXP ZR, SEXP weightsR);
RcppExport SEXP posKernel(SEXP posR, SEXP widthR);

#endif
