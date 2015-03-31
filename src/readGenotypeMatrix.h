//
// File readGenotypeMatrix.h defining prototypes for readGenotypeMatrix.cpp
//

#ifndef _readGenotypeMatrix_H_

#define _readGenotypeMatrix_H_

RcppExport SEXP readGenotypeMatrix(SEXP ext, SEXP seqnamesR, SEXP startR,
				   SEXP endR, SEXP subsetR,
				   SEXP noIndelsR, SEXP onlyPassR,
				   SEXP NAlimitR, SEXP MAFlimitR,
				   SEXP NAactionR, SEXP MAFactionR,
                                   SEXP sexR);

#endif
