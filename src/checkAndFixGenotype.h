//
// File checkAndFixGenotype.h defining prototypes for checkAndFixGenotype.cpp
//

#ifndef _checkAndFixGenotype_H_

#define _checkAndFixGenotype_H_

RcppExport SEXP checkAndFixGenotype(SEXP iR, SEXP pR, SEXP DimR, SEXP xR,
				    SEXP ploidyR,
				    SEXP NAlimitR, SEXP MAFlimitR,
                                    SEXP NAactionR, SEXP MAFactionR);
RcppExport SEXP checkAndFixGenotypeChar(SEXP ZR,
					SEXP NAlimitR, SEXP MAFlimitR,
					SEXP NAactionR, SEXP MAFactionR,
                                        SEXP sexR);

#endif
