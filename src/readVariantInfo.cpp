#include <R.h>
#include <Rcpp.h>

#include <vector>
#include <string>

#include <htslib/tbx.h>
#include <htslib/bgzf.h>

using namespace Rcpp;
using namespace std;

#include "errCodes.h"
#include "readVariantInfo.h"

#define TYPE_INDEL          1
#define TYPE_MULTIPLE       2
#define TYPE_TRANSITION     3
#define TYPE_TRANSVERSION   4
#define TYPE_UNKNOWN        5


/*-------------------------------------------------------------------------*/
/* Some definitions copied from Rsamtools                                  */
/*-------------------------------------------------------------------------*/

typedef struct {
    htsFile *file;
    tbx_t *index;
    hts_itr_t *iter;
} _TABIX_FILE;

#define TABIXFILE(b) ((_TABIX_FILE *) R_ExternalPtrAddr(b))

static int64_t _tbx_tell(htsFile *file)
{
    BGZF *fp;

    if (!file->is_bgzf)
        Rf_error("[internal] hmm.. this doesn't look like a tabix file, sorry");
    fp = file->fp.bgzf;
    return bgzf_tell(fp);
}
static void _tbx_seek(htsFile *file, int64_t offset)
{
    BGZF *fp;

    if (!file->is_bgzf)
        Rf_error("[internal] hmm.. this doesn't look like a tabix file, sorry");
    fp = file->fp.bgzf;
    if (bgzf_seek(fp, offset, SEEK_SET) < 0)
        Rf_error("[internal] bgzf_seek() failed");
    return;
}
static const char *_tbx_read_line(htsFile *file, int *len)
{
    BGZF *fp;
    static kstring_t ksbuf = {0, 0, NULL};

    if (!file->is_bgzf)
        Rf_error("[internal] hmm.. this doesn't look like a tabix file, sorry");
    fp = file->fp.bgzf;
    if (bgzf_getline(fp, '\n', &ksbuf) < 0)
        return NULL;
    *len = ksbuf.l;
    return ksbuf.s;
}
static const char *_tbx_read_next_rec(htsFile *file, tbx_t *index,
                                      hts_itr_t *iter, int *len)
{
    static kstring_t ksbuf = {0, 0, NULL};

    if (tbx_itr_next(file, index, iter, &ksbuf) < 0)
        return NULL;
    *len = ksbuf.l;
    return ksbuf.s;
}
static void _skip_header_lines(htsFile *file, const tbx_conf_t *conf)
{
    const char *line;
    int linelen;

    int64_t curr_off = _tbx_tell(file);
    while ((line = _tbx_read_line(file, &linelen)) != NULL) {
        if (line[0] != conf->meta_char)
            break;
        curr_off = _tbx_tell(file);
    }
    _tbx_seek(file, curr_off);
    return;
}

/*- end of definitions and functions copied from Rsamtools ----------------*/


/*-------------------------------------------------------------------------*/
/* function that determines type of variant from ref and alt descriptors   */
/*-------------------------------------------------------------------------*/

int determineType(string &ref, string &alt)
{

    if (ref.length() == 1 && alt.length() == 1)
    {
	if (ref == "A")
	{
	    if (alt == "G")
		return TYPE_TRANSITION;
	    else
		return TYPE_TRANSVERSION;		    
	}
	else if (ref == "T")
	{
	    if (alt == "C")
		return TYPE_TRANSITION;
	    else
		return TYPE_TRANSVERSION;		    
	}
	else if (ref == "G")
	{
	    if (alt == "A")
		return TYPE_TRANSITION;
	    else
		return TYPE_TRANSVERSION;		    
	}
	else // if (ref == "C")
	{
	    if (alt == "T")
		return TYPE_TRANSITION;
	    else
		return TYPE_TRANSVERSION;		    
	}
    }
    else if (ref.length() > 1)
	return TYPE_INDEL;
    else
    {
	int maxLen = 0, noComma = 0, len = 0;

	// go through ALT field
	for (string::iterator it = alt.begin(); it != alt.end(); ++it)
	{
	    if (*it == ',') // switch alternative allele
	    {
		if (len > maxLen)
		    maxLen = len;

		len = 0;
		noComma++;
	    }
	    else
		len++;
	}
	    
	if (len > maxLen)
	    maxLen = len;

	if (maxLen > 1)
	    return TYPE_INDEL;
	else if (maxLen == 1 && noComma >= 1)
	    return TYPE_MULTIPLE;
	else
	    return TYPE_UNKNOWN;
    }
}


/*-------------------------------------------------------------------------*/
/* function that reads record by record of a region and adds the data to   */
/*     vectors of descriptors                                              */
/*-------------------------------------------------------------------------*/

int tabixToInfo(htsFile *file, tbx_t *index, hts_itr_t *iter,
		int subsetLen, int *subset, int sexLen, int *sex,
		bool noIndels, bool onlyPass,
		double NAlimit, double MAFlimit, int NAaction, int MAFaction,
		bool refAlt, bool omitZeroMAF,
		vector<string> &names, vector<string> &chrom,
		vector<int> &pos, vector<int> &type,
		vector<string> &refV, vector<string> &altV,
		vector<double> &MAF, int &n, int &m)
{
    int linelen;
    const char *line;
    const tbx_conf_t conf = index->conf;
    double colMAF, colNA;
    vector<int> entryBuffer, ploidy;

    // read from tabix stream record by record
    while ((line = _tbx_read_next_rec(file, index, iter, &linelen)) != NULL) {
        if (line[0] == conf.meta_char)
            continue;

	int spos;
	string seqname, id, ref, alt, qual, filter, info, format;
	istringstream lineBuffer(line);

	// read fixed fields
	lineBuffer >> seqname >> spos >> id >> ref >> alt >> qual >> filter
	           >> info >> format;

	if (lineBuffer.fail())
	    return ERR_ROWSYNTAX;

	// check for input filter and GT format field
	if ((onlyPass && filter != "PASS") || format.substr(0, 2) != "GT")
	    continue;
	else if (noIndels) // check for indels
	{
	    if (ref.length() != 1) // REF longer than 1 => move to next
		continue;

	    int current = 0, maxi = 0;

	    for (string::iterator it = alt.begin(); it != alt.end(); ++it)
	    {
		if (*it == ',')
		    current = 0;
		else
		    current++;

		if (current > maxi)
		    maxi = current;
	    }
	    
	    if (maxi > 1) // at least one ALT longer than 1 => move to next
		continue;
	}

	int col = 0, vcol = 0, colSum = 0, noNAs = 0, totalChrs = 0;

	// read genotype columns
	while(lineBuffer.good())
	{
	    string entry;

	    lineBuffer >> entry; // next column

	    if (subsetLen > 0)
	    {
		if (vcol >= subsetLen)
		    break;
		else if (subset[vcol] == 0)
		{
		    vcol++;
		    continue;
		}
	    }

	    int current = 0, total = 0, currentPloidy = 0;
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
			return ERR_GTINVALID;

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
		else if (*it >= '1' && *it <= '9')
		    current = 1;
		else if (*it != '0') // invalid character
		    return ERR_GTINVALID;
	    }

	    if (waitingForGT) // GT incomplete
		return ERR_GTINVALID;

	    total += current;
	    currentPloidy++;
	    totalChrs += currentPloidy;
	    colSum += (total % 256);
	    noNAs += (total / 256);

	    if (total == 1 && currentPloidy == 2 && sexLen > 0 &&
		sex[col] > 1)
		colSum++;

	    col++;
	    vcol++;
	}

	if (!col) // nothing read => move on to next
	    continue;

	if (m == 0)
	    n = col;
	else if (col != n) // check against no. of columns
	    return ERR_NOMISMATCH;

	if (noNAs && NAaction == 4)
	    return ERR_MISSING;
	else if (!noNAs || NAaction != 2)
	{
	    colMAF = colSum * 1. / totalChrs;
	    colNA = noNAs * 1. / totalChrs;

	    if (colMAF > 0.5 && MAFaction == 4)
		return ERR_MAF50;
	    else if (colNA <= NAlimit)
	    {
		if (colMAF > 0.5 && MAFaction == 1) // invert column
		    colMAF = 1. - colMAF;

		if (colMAF <= MAFlimit &&
		    (colMAF <= 0.5 || MAFaction == 3) && 
		    (colMAF > 0 || !omitZeroMAF))
		{
		    chrom.push_back(seqname);
		    pos.push_back(spos);
		    MAF.push_back(colMAF);

		    if (id == ".")
		    {
			ostringstream outId;
			outId << seqname << ":" << spos;
			names.push_back(outId.str());
		    }
		    else
			names.push_back(id);

		    if (refAlt)
		    {
			refV.push_back(ref);
			altV.push_back(alt);
		    }

		    type.push_back(determineType(ref, alt));
		    
		    m++;
		}
	    }
	}
    }

    return ERR_OK;
}


/*-------------------------------------------------------------------------*/
/* function that reads a VCF file and stores the data into a sparse matrix */
/*     (multiple regions are merged into one single matrix)                */
/*-------------------------------------------------------------------------*/

RcppExport SEXP readVariantInfo(SEXP ext, SEXP seqnamesR, SEXP startR,
				SEXP endR, SEXP subsetR,
				SEXP noIndelsR, SEXP onlyPassR,
				SEXP NAlimitR, SEXP MAFlimitR,
				SEXP NAactionR, SEXP MAFactionR,
				SEXP refAltR, SEXP omitZeroMAFR, SEXP sexR)
{
    CharacterVector seqnames(seqnamesR);
    IntegerVector start(startR), end(endR), sex(sexR);
    LogicalVector subset(subsetR);
    int nspc = seqnames.length(), n = 0, m = 0;
    bool noIndels = as<bool>(noIndelsR), onlyPass = as<bool>(onlyPassR);
    double NAlimit = as<double>(NAlimitR);
    double MAFlimit = as<double>(MAFlimitR);
    int NAaction = as<int>(NAactionR);
    int MAFaction = as<int>(MAFactionR);
    bool omitZeroMAF = as<bool>(omitZeroMAFR);
    bool refAlt = as<bool>(refAltR);
    _TABIX_FILE *tfile = TABIXFILE(ext);
    htsFile *file = tfile->file;
    tbx_t *index = tfile->index;
    vector<string> chrom, names, ref, alt;
    vector<int> pos, type;
    vector<double> MAF;

    if (nspc == 0) // no region defined => read entire VCF file
    {
        hts_itr_t *iter = tfile->iter;
        if (iter == NULL) {
            _skip_header_lines(file, &(index->conf));
            /* Create iterator that will iterate over all records. */
            iter = tbx_itr_queryi(index, HTS_IDX_REST, 0, 0);
            if (iter == NULL)
                Rf_error("[internal] failed to create tabix iterator");
            tfile->iter = iter;
        }
	// read VCF into sparse matrix
        int ret = tabixToInfo(file, index, iter,
                              subset.length(),
			      (subset.length() > 0 ? subset.begin() : NULL),
			      sex.length(),
			      (sex.length() > 0 ? sex.begin() : NULL),
			      noIndels, onlyPass, 
			      NAlimit, MAFlimit, NAaction, MAFaction, refAlt,
			      omitZeroMAF, names, chrom, pos, type, ref, alt,
			      MAF, n, m);

	if (ret != ERR_OK)
	    return ERR_RETURN(ret);
    }
    else // region(s) specified
    {
        for (int ispc = 0; ispc < nspc; ++ispc) // iterate over regions
	{
            int ibeg, iend, tid;
            hts_itr_t *iter;

            ibeg = start[ispc] == 0 ? 0 : start[ispc] - 1;
            iend = end[ispc];

            tid = tbx_name2id(index, seqnames[ispc]);
            if (tid < 0)
                return ERR_RETURN(ERR_CHRMISSING);

            iter = tbx_itr_queryi(index, tid, ibeg, iend);

	    // read VCF into matrix (columns accumulate over multiple regions)
	    int ret = tabixToInfo(file, index, iter,
                                  subset.length(),
				  (subset.length() > 0 ? subset.begin() : NULL),
				  sex.length(),
				  (sex.length() > 0 ? sex.begin() : NULL),
				  noIndels, onlyPass, 
				  NAlimit, MAFlimit, NAaction, MAFaction,
				  refAlt, omitZeroMAF,
				  names, chrom, pos, type, ref, alt,
				  MAF, n, m);

            tbx_itr_destroy(iter);

	    if (ret != ERR_OK)
		return ERR_RETURN(ret);
         }
    }

    // return NULL if no data have been read
    if (!m)
	return R_NilValue;

    // return list with output
    List ret;
    
    ret["names"] = wrap(names);
    ret["seqnames"] = wrap(chrom);
    ret["pos"] = wrap(pos);
    ret["MAF"] = wrap(MAF);

    if (refAlt)
    {
	ret["ref"] = wrap(ref);
	ret["alt"] = wrap(alt);
    }

    CharacterVector typeLevels =
	CharacterVector::create("INDEL", "MULTIPLE", "TRANSITION",
				"TRANSVERSION", "UNKNOWN");
    IntegerVector typeFactor = wrap(type);
    typeFactor.attr("levels") = typeLevels;
    typeFactor.attr("class") = "factor";

    ret["type"] = typeFactor;
    
    return ret;
}
