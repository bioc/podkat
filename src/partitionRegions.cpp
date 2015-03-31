#include <Rcpp.h>

#include <vector>
#include <string>

using namespace Rcpp;
using namespace std;

#include "partitionRegions.h"

/*-------------------------------------------------------------------------*/
/* function for generating a partition of overlapping windows for a given  */
/*     list of windows                                                     */
/*-------------------------------------------------------------------------*/

RcppExport SEXP partitionRegions(SEXP seqNamesR, SEXP startR, SEXP endR,
                                 SEXP widthR, SEXP overlapR)
{
    CharacterVector seqNames(seqNamesR);
    IntegerVector start(startR), end(endR);
    int width = as<int>(widthR), n = start.length();
    double overlap = as<double>(overlapR);
    int shift = (int)floor(width - width * overlap);
    vector<string> outSeqNames;
    vector<int> outStart, outEnd;

    for (int i = 0; i < n; i++)
    {
	int s = start[i], e = end[i], l = e - s + 1;
	string chr(seqNames[i]);

	if (l <= width)
	{
	    outSeqNames.push_back(chr);
	    outStart.push_back(s);
	    outEnd.push_back(e);
	}
	else
	{
	    int m = ceil((l - width) / (1. * shift));
	    int overhang = s + m * shift + width - 1 - e;
	    int offset = -(int)round(overhang / 2.);

	    outSeqNames.push_back(chr);
	    outStart.push_back(s);
	    outEnd.push_back(s + width + offset - 1);

	    for (int j = 1; j < m; j++)
	    {
		outSeqNames.push_back(chr);
		outStart.push_back(s + offset + j * shift);
		outEnd.push_back(s + offset + j * shift + width - 1);
	    }

	    outSeqNames.push_back(chr);
	    outStart.push_back(s + offset + m * shift);
	    outEnd.push_back(e);
	}
    }

    List ret;

    ret["seqnames"] = wrap(outSeqNames);
    ret["start"] = wrap(outStart);
    ret["end"] = wrap(outEnd);

    return ret;
}
