#include <R_ext/Rdynload.h>
#include <Rcpp.h>

using namespace Rcpp;

#include "pValues.h"
#include "kernels.h"
#include "bernoulliExact.h"
#include "cumMax.h"
#include "doubleMale.h"
#include "checkAndFixGenotype.h"
#include "partitionRegions.h"
#include "readGenotypeMatrix.h"
#include "readVariantInfo.h"


static const R_CallMethodDef callMethods[] = 
{
    /* pValues.cpp */
    {"davies", (DL_FUNC) &davies, 3},
    {"liu", (DL_FUNC) &liu, 2},
    {"liuMod", (DL_FUNC) &liuMod, 2},
    /* kernels.cpp */
    {"localSimKernel", (DL_FUNC) &localSimKernel, 1},
    {"localSimKernelWeighted", (DL_FUNC) &localSimKernelWeighted, 2},
    {"posKernel", (DL_FUNC) &posKernel, 2},
    /* bernoulliExact.cpp */
    {"computeExactBernoulliPvalue", (DL_FUNC) &computeExactBernoulliPvalue, 4},
    /* cumMax.cpp */
    {"cumMax", (DL_FUNC) &cumMax, 2},
    /* doubleMale.cpp */
    {"doubleMale", (DL_FUNC) &doubleMale, 3},
    /* checkAndFixGenotype.cpp */
    {"checkAndFixGenotype", (DL_FUNC) &checkAndFixGenotype, 9},
    {"checkAndFixGenotypeChar", (DL_FUNC) &checkAndFixGenotypeChar, 6},
    /* partitionRegions.cpp */
    {"partitionRegions", (DL_FUNC) &partitionRegions, 5},
    /* readGenotypeMatrix.cpp */
    {"readGenotypeMatrix", (DL_FUNC) &readGenotypeMatrix, 12},
    /* readVariantInfo.cpp */
    {"readVariantInfo", (DL_FUNC) &readVariantInfo, 14},
    {NULL, NULL, 0}
};

extern "C" 
{

    void R_init_podkat(DllInfo *info)
    {
        /* register routines, allocate resources */
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
    }

    void R_unload_podkat(DllInfo *info)
    {
        /* release resources */
    }
}
