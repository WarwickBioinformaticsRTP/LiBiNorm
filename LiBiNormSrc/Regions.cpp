#include <vector>
#include "Regions.h"
#include "parser.h"

using namespace std;
using namespace BamTools;

//	Combines a region with an existing set of regions
void regionList::combineRegion(const region & r1)
{
	bool combined = false;

	for (auto i = data.begin();i != data.end();i++)
	{
		//	Does the new region start before the end of the existinng region, and is it on the same strand?
		if ((r1.start <= i->second.end) && (r1.strand == i->second.strand))
		{
			//	If it starts before the existing start and ends after it then there is an overlap 
			//	and we are going to have to replace the existing region as regions are indexed by
			//  the start
			if ((r1.start <= i->second.start) && (r1.end >= i->second.start)) 
			{
				// The new region, based on the existing one
				region r = i->second;
				// but with an earlier start
				r.start = r1.start;
				//	and possibly an earlier end.
				if (r1.end >= i->second.end)
					r.end = r1.end;

				//	replace the existing entry with the new one, taking care not to saw off the
				//	branch you are sitting on and 
				data.erase(i);
				i = data.emplace(r1.start,r).first;
				combined = true;
				break;	
			}
			//	The new region just estends the end of the existing region
			else if (r1.end >= i->second.end)
			{
				i->second.end = r1.end;
				combined = true;
				break;	
			}
		}
	}
	if (!combined)
		data.emplace(r1.start,r1);
}

//	Combine this region list with another (they are on the same chromosome
void regionList::combine(const regionList & rl)
{
	for (auto r : rl.data)
		combineRegion(r.second);
}


regionList::regionList(const readData & read)
{
	// initialize alignment end to starting position

	size_t start = read.position;
	size_t end = start;

	// iterate over cigar operations
	vector<CigarOp>::const_iterator cigarIter = read.cigar.begin();
	vector<CigarOp>::const_iterator cigarEnd  = read.cigar.end();
	for ( ; cigarIter != cigarEnd; ++cigarIter) {
		const CigarOp& op = (*cigarIter);

		switch ( op.Type ) {

			// increase end position on CIGAR chars [DMXN=]
			case Constants::BAM_CIGAR_DEL_CHAR      :
			case Constants::BAM_CIGAR_MATCH_CHAR    :
			case Constants::BAM_CIGAR_MISMATCH_CHAR :
			case Constants::BAM_CIGAR_SEQMATCH_CHAR :
				end += op.Length;
				break;

			case Constants::BAM_CIGAR_INS_CHAR :
				break;

			case Constants::BAM_CIGAR_REFSKIP_CHAR  :
				{
					combineRegion(region(start,end-1,read.strand));
					start = (end + op.Length);
					end = start;
					break;
				}

			default :
				break;
		}
	}

	combineRegion(region(start,end-1,read.strand));
}

//	Parse a text string and convert it into a cigar value.  Used when retreiving entries from cache files
void parserInternal::parseval(const char *& start,Cigar & co,size_t & len)
{
	int i = 0;
	while (i < len)
	{
		char c = start[i++];
		int val = 0;
		while ((start[i] >= '0') && (start[i] <= '9') && (i < len))
			val = (val *10)+ (start[i++]-'0');
		co.emplace_back(CigarOp(c,val));
	}
}

//	Prints out cigar information
bool printVal(outputDataFile * f,const Cigar & cigar)
{
	for (auto & i: cigar)
		fprintf(f->fout,"%c%i",i.Type,i.Length);
	return true;
};

//	Prints the value of a readData instance.  Note that printing read.cigar will have the effect of calling the 
//	printVal associated with the cigar above
bool printVal(outputDataFile * f,const readData & read)
{
	f->printStart(printZero(read.refId),read.position,read.strand,read.cigar);
	return true;
};


