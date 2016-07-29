#include <vector>
#include "Regions.h"
#include "parser.h"

using namespace std;
using namespace BamTools;

void regionList::combineRegion(const region & r1)
{
	bool combined = false;

	for (auto i = data.begin();i != data.end();i++)
	{
		if ((r1.start <= i->second.end) && (r1.strand == i->second.strand))
		{
			if ((r1.start <= i->second.start) && (r1.end >= i->second.start)) 
			{
				region r = i->second;
				r.start = r1.start;
				if (r1.end >= i->second.end)
					r.end = r1.end;

				if (i == data.begin())
				{
					data.erase(i);
					data.emplace(r1.start,r);
					i = data.begin();
				}
				else
				{
					auto j = next(i,-1);
					data.erase(i);
					data.emplace(r1.start,r);
					i = next(j,1);
				}
				combined = true;
				break;	
			}
			if (r1.end >= i->second.end)
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

void regionList::combine(const regionList & rl)
{
	for (const auto r : rl.data)
		combineRegion(r.second);
}


void regionList::add(int start,int end,char strand)
{
	data.emplace(start,region(start,end,strand));
}

void regionLists::GetRegions(const BamAlignment & ba) 
{ 	
	data[ba.RefID].GetRegions(ba);
}

void regionLists::combine(const regionLists & rl)
{
	for (auto i : rl.data)
		data[i.first].combine(i.second);
}

/*
regionList::regionList(const cacheEntry & read): regionList::regionList(read.start,
{

}
*/


void regionList::GetRegions(const BamAlignment & ba) {

//	bool revStrand = (ba.IsReverseStrand() == ba.IsFirstMate());

//	auto & CigarData = ba.CigarData;
//	size_t start = ba.Position +1;

	cacheEntry  cr(ba);
//	regionList rl(cr);

	combine(regionList(cr));
}

regionList::regionList(size_t start,char strand,const std::vector<BamTools::CigarOp> & CigarData)
{
	// initialize alignment end to starting position

	size_t end = start;

	// iterate over cigar operations
	vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
	vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
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
					combineRegion(region(start,end-1,strand));
					start = (end + op.Length);
					end = start;
					break;
				}

			default :
				break;
		}
	}

	combineRegion(region(start,end-1,strand));
}

/*
regionLists::regionLists(const string & line)
{
	parseTsv(line,name,data);
}
*/
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
//	parseTsv(start,rl.data);
}

/*
void parserInternal::parseval(const char *& start,regionList & rl,size_t & len)
{
	parseTsv(start,rl.data);
}
void parserInternal::parseval(const char *& start,region & r,size_t & len)
{
	parseTsv(start,r.start,r.end,r.strand);
}
*/

