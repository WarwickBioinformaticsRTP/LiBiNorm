#include <vector>
#include "Regions.h"

using namespace std;
using namespace BamTools;

bool regionList::combineRegion(size_t start,size_t finish,bool revStrand)
{
	bool combined = false;

	for (iterator i = begin();i != end();i++)
	{
		if (start <= i->second.end)
		{
			if ((start <= i->second.start) && (finish >= i->second.start)) 
			{
				region r = i->second;
				r.start = start;
				if (finish >= i->second.end)
					r.end = finish;

				if (i == begin())
				{
					erase(i);
					emplace(start,r);
					i = begin();
				}
				else
				{
					iterator j = next(i,-1);
					erase(i);
					emplace(start,r);
					i = next(j,1);
				}
				combined = true;
				break;	
			}
			if (finish >= i->second.end)
			{
				i->second.end = finish;
				combined = true;
				break;	
			}
		}
	}
	if (!combined)
	{
		emplace(start,region(start,finish,revStrand));
		return false;
	}
	return true;
}

void regionList::add(size_t start,size_t end,bool revStrand)
{
	emplace(start,region(start,end,revStrand));
}

void regionLists::GetRegions(const BamAlignment & ba) { 	This[ba.RefID].GetRegions(ba);
}


void regionList::GetRegions(const BamAlignment & ba) {

	//	If we already have some regions then we need to combine them
	bool combine = size();
//	standardPair = true;

	_DBG(
		bool firstMate = ba.IsFirstMate();
	bool revSt = ba.IsReverseStrand();
	bool mateMapped = ba.IsMateMapped();

	bool paired = ba.IsProperPair();
	)

/*	bool revStrand = ba.IsMateMapped()?
		(ba.IsProperPair()?(ba.IsReverseStrand() == ba.IsFirstMate()):!ba.IsReverseStrand()):
		ba.IsReverseStrand();
*/
	bool revStrand = (ba.IsReverseStrand() == ba.IsFirstMate());

	if (!ba.IsMateMapped())
	{
		revStrand = (ba.IsReverseStrand() == ba.IsFirstMate());
//		revStrand = ba.IsReverseStrand();
	}
	else if (!ba.IsProperPair())
	{
		revStrand = (ba.IsReverseStrand() == ba.IsFirstMate());
//		standardPair = false;
//		revStrand = !ba.IsReverseStrand();
	}
	
	//	Dont combine if the two reads are on reverse strands
	if (combine && (revStrand != (begin()->second.strand == '-')))
		combine = false;

	auto & CigarData = ba.CigarData;

	// initialize alignment end to starting position

	size_t start = ba.Position +1;
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
					if (combine)
						combineRegion(start,end-1,revStrand);
					else
						add(start,end-1,revStrand);
					end = start = (end + op.Length);
					start++;		//Not convinced that the increement should be here, but is required for
									//compatibility with htseq-count.
					break;
				}

			default :
				break;
		}
	}
	if (combine)
		combineRegion(start,end-1,revStrand);
	else
		add(start,end-1,revStrand);
}

