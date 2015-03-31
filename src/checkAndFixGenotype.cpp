#include <limits.h>
#include <float.h>
#include <Rcpp.h>

#include <vector>

using namespace Rcpp;
using namespace std;

#include "errCodes.h"
#include "checkAndFixGenotype.h"


RcppExport SEXP checkAndFixGenotype(SEXP iR, SEXP pR, SEXP DimR, SEXP xR,
				    SEXP ploidyR,
				    SEXP NAlimitR, SEXP MAFlimitR,
                                    SEXP NAactionR, SEXP MAFactionR)
{
    IntegerVector i(iR);
    IntegerVector p(pR);
    IntegerVector Dim(DimR);
    NumericVector x(xR);
    double NAlimit = as<double>(NAlimitR);
    double MAFlimit = as<double>(MAFlimitR);
    int NAaction = as<int>(NAactionR);
    int MAFaction = as<int>(MAFactionR);
    int ploidy = as<int>(ploidyR);
    int noR = Dim[0], noC = Dim[1];
    int J, Index, noNAs, noEntries = 0, newNoCols = 0;
    double colSum, colNA, colMAF;
    vector<double> newX, MAF;
    vector<int> newP, newI;
    LogicalVector colOmitted(noC);

    newP.push_back(0);

    for (J = 0; J < noC; J++)
    {
	colSum = 0;
	noNAs = 0;
	colOmitted[J] = 1;

	for (Index = p[J]; Index < p[J + 1]; Index++)
	{
	    if (R_IsNA(x[Index]))
		noNAs++;
	    else if (x[Index] > 0 && x[Index] <= ploidy)
	    {
		colSum += x[Index];
		noEntries++;
	    }
	    else if (x[Index] != 0)
		return ERR_RETURN(ERR_INVALIDMAT);
	}

	if (noNAs && NAaction == 4)
	    return ERR_RETURN(ERR_MISSINGMAT);
	else if (noEntries && (!noNAs || NAaction != 2))
	{
	    colMAF = colSum / (ploidy * noR);
	    colNA = noNAs * 1. / noR;

	    if (colMAF > 0.5 && MAFaction == 4)
		return ERR_RETURN(ERR_MAF50);
	    else if (colNA <= NAlimit)
	    {
		int EntriesAdded = 0;

		if (colMAF > 0.5 && MAFaction == 1)
		{
		    NumericVector aux(noR);

		    for (Index = p[J]; Index < p[J + 1]; Index++)
			if (R_IsNA(x[Index]))
			    aux[i[Index]] = ploidy;
			else
			    aux[i[Index]] = x[Index];

		    colMAF = 0;

		    for (int I = 0; I < noR; I++)
		    {
			aux[I] = ploidy - aux[I];
			colMAF += aux[I];
		    }

		    colMAF /= (ploidy * noR);

		    if (colMAF <= MAFlimit && colMAF > 0)
		    {
			for (int I = 0; I < noR; I++)
			{
			    if (aux[I] > 0)
			    {
				newI.push_back(I);
				newX.push_back(aux[I]);
				EntriesAdded++;
			    }
			}
		    }
		}
		else if (colMAF <= MAFlimit &&
			 (colMAF <= 0.5 || MAFaction == 3))
		{
		    for (Index = p[J]; Index < p[J + 1]; Index++)
		    {
			if (!R_IsNA(x[Index]))
			{
			    newI.push_back(i[Index]);
			    newX.push_back(x[Index]);
			    EntriesAdded++;
			}
		    }
		}
		
		if (EntriesAdded)
		{
		    newP.push_back(newP[newP.size() - 1] + EntriesAdded);
		    MAF.push_back(colMAF);
		    newNoCols++;
		    colOmitted[J] = 0;
		}
	    }
	}
    }

    List ret;

    ret["i"] = wrap(newI);
    ret["p"] = wrap(newP);
    ret["x"] = wrap(newX);
    ret["Dim"] = IntegerVector::create(noR, newNoCols);
    ret["omitted"] = colOmitted;
    ret["MAF"] = wrap(MAF);

    return(ret);
}


RcppExport SEXP checkAndFixGenotypeChar(SEXP ZR,
					SEXP NAlimitR, SEXP MAFlimitR,
					SEXP NAactionR, SEXP MAFactionR,
                                        SEXP sexR)
{
    CharacterMatrix Z(ZR);
    IntegerVector sex(sexR);
    double NAlimit = as<double>(NAlimitR);
    double MAFlimit = as<double>(MAFlimitR);
    int NAaction = as<int>(NAactionR);
    int MAFaction = as<int>(MAFactionR);
    int n = Z.nrow(), m = Z.ncol();
    int I, J, noNAs, colSum, newNoCols = 0;
    double colMAF, colNA;
    vector<double> x, MAF;
    vector<int> i, p;
    IntegerVector entryBuffer(n), ploidy(n);
    LogicalVector colOmitted(m);

    p.push_back(0); // p field of dgCMatrix has to start with 0

    for (J = 0; J < m; J++) // proceed by column
    {
	int totalChrs = 0;
	colOmitted[J] = 1;

	for (I = 0; I < n; I++)
	{
	    string entry = (char *)Z(I, J);
	    int current = 0, total = 0, currentPloidy = 1;
	    bool waitingForGT = true; // what character to expect

	    // go through GT field of current column
	    for (string::iterator it = entry.begin(); it != entry.end(); ++it)
	    {
		if (waitingForGT)
		{
		    if (*it >= '1' && *it <= '9') // non-zero field found
			current = 1;
		    else if (*it == '.') // missing value
			current = 256;
		    else if (*it != '0') // invalid character
			return ERR_RETURN(ERR_GTINVALID);

		    waitingForGT = false;
		}
		else if (*it == '/' || *it == '|') // switch haplotype
		{
		    total += current;
		    current = 0;
		    currentPloidy++;

		    waitingForGT = true;
		}
		else if (*it == ':') // end of GT field
		{
		    waitingForGT = false;
		    break;
		}
		else // invalid character
		    return ERR_RETURN(ERR_GTINVALID);
	    }

	    if (waitingForGT) // GT incomplete
		return ERR_RETURN(ERR_GTINVALID);
	    
	    total += current;
	    entryBuffer[I] = total;
	    ploidy[I] = currentPloidy;

	    if (entryBuffer[I] == 1 && ploidy[I] == 2 && sex.length() > 0 &&
		sex[I] > 1)
		entryBuffer[I] = 2;

	    totalChrs += currentPloidy;
	}

	colSum = 0;
	noNAs = 0;

	// determine number of NAs and non-zeros
	for (I = 0; I < n; I++)
	{
	    if (entryBuffer[I] >= 256)
	    {
		colSum += (entryBuffer[I] % 256);
		noNAs += (entryBuffer[I] / 256);
	    }
	    else
		colSum += entryBuffer[I];
	}

	if (noNAs && NAaction == 4)
	    return ERR_RETURN(ERR_MISSINGMAT);
	else if (colSum && (!noNAs || NAaction != 2))
	{
	    colMAF = colSum * 1. / totalChrs;
	    colNA = noNAs * 1. / totalChrs;

	    if (colMAF > 0.5 && MAFaction == 4)
		return ERR_RETURN(ERR_MAF50);
	    else if (colNA <= NAlimit)
	    {
		int EntriesAdded = 0;

		if (colMAF > 0.5 && MAFaction == 1) // invert column
		{
		    if (noNAs + colSum < totalChrs)
		    {
			NumericVector aux(n);
			
			colMAF = 0;
			
			for (I = 0; I < n; I++)
			{
			    aux[I] = ploidy[I] - (entryBuffer[I] % 256)
				               - (entryBuffer[I] / 256);

			    colMAF += aux[I];
			}

			colMAF /= totalChrs;

			if (colMAF <= MAFlimit && colMAF > 0)
			{
			    for (I = 0; I < n; I++)
			    {
				if (aux[I] > 0)
				{
				    x.push_back(aux[I]);
				    i.push_back(I);
				    EntriesAdded++;
				}
			    }
			}
		    }
		}
		else if (colMAF <= MAFlimit &&
			 (colMAF <= 0.5 || MAFaction == 3))
		{
		    for (I = 0; I < n; I++)
		    {
			int val = entryBuffer[I] % 256;

			if (val)
			{
			    x.push_back((double)val);
			    i.push_back(I);
			    EntriesAdded++;
			}
		    }
		}
		else
		    EntriesAdded = 0;
		
		if (EntriesAdded)
		{
		    // fill up output at end of row
		    p.push_back(p[p.size() - 1] + EntriesAdded);
		    MAF.push_back(colMAF);
		    newNoCols++;
		    colOmitted[J] = 0;
		}
	    }
	}
    }

    List ret;

    ret["i"] = wrap(i);
    ret["p"] = wrap(p);
    ret["x"] = wrap(x);
    ret["Dim"] = IntegerVector::create(n, newNoCols);
    ret["omitted"] = colOmitted;
    ret["MAF"] = wrap(MAF);

    return(ret);
}
