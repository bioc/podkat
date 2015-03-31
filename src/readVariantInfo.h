//
// File readVariantInfo.h defining prototypes for readVariantInfo.cpp
//

#ifndef _readVariantInfo_H_

#define _readVariantInfo_H_

RcppExport SEXP readVariantInfo(SEXP ext, SEXP seqnamesR, SEXP startR,
				SEXP endR, SEXP subsetR,
				SEXP noIndelsR, SEXP onlyPassR,
				SEXP NAlimitR, SEXP MAFlimitR,
				SEXP NAactionR, SEXP MAFactionR,
				SEXP refAltR, SEXP omitZeroMAFR, SEXP sexR);

#endif
