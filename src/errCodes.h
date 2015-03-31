//
// File errCodes.h defining error codes and messages
//

#ifndef _errCodes_H_

#define _errCodes_H_

/* error handling */

#define ERR_OK         0
#define ERR_READLINE   1
#define ERR_ROWSYNTAX  2
#define ERR_GTINVALID  3
#define ERR_INDEXFAIL  4
#define ERR_NOMISMATCH 5
#define ERR_MISSING    6
#define ERR_MAF50      7
#define ERR_CHRMISSING 8
#define ERR_MISSINGMAT 9
#define ERR_INVALIDMAT 10

static const char *errCodes[] =
{
    "",
    "read line failed; corrupt or invalid file?",
    "row in VCF file does not conform to standard format",
    "invalid character in GT field",
    "failed to load tabix index",
    "inconsistent numbers of samples in VCF file",
    "missing values in VCF file",
    "MAF > 0.5 in genotype matrix",
    "requested sequence name not present in tabix index",
    "missing values in genotype matrix",
    "invalid values in genotype matrix",
};

#define ERR_RETURN(CODE) CharacterVector::create(errCodes[CODE])

#endif
