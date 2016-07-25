#include "GtfFileEx.h"

using namespace std;

void gtfRegion::checkOverlap(const region & segment,vector<gtfOverlap> & overlaps) const
{
	//	strict fit
	if ((segment.start >= start) && (segment.end <= finish))
	{
		overlaps.emplace_back(segment.start,segment.end,true,name,type);
	}
	//	non-strict fits
	else if ((segment.start <= finish) && (segment.start >= start)) 
	{
		overlaps.emplace_back(segment.start,finish,false,name,type);
	}
	else if ((segment.end <= finish) && (segment.end >= start))
	{
		overlaps.emplace_back(start,segment.end,false,name,type);
	}
	else if ((segment.start <= start) && (segment.end >= finish))
	{
		overlaps.emplace_back(start,finish,false,name,type);
	}
}



void gtfFileEx::index(geneCountsClass & geneCounts,const std::string & type)
{

	for (auto & chrom : entryMap)
	{
		//	In each chromosome
		chromosomeGtfData & thisChromData = genomeGtfData[chrom.first];
		for (auto i = chrom.second.begin(); i != chrom.second.end();i++)
		{
			size_t finish = i->second.finish;
			setEx<string> type(i->second.type);
			for (auto j = next(i,1);(j != chrom.second.end()) && (j->first <= (finish + 1));)
			{
				auto k = j++;
				if (k != chrom.second.end())
				{
					if ((i->second.type == k->second.type) && (i->second.tags[0].val == k->second.tags[0].val))
					{
						if (k->second.finish > finish)
						{
							finish = k->second.finish;
							type.add(k->second.type);
//							k->second.valid = false;
							chrom.second.erase(k);
						}
						else if (k->second.finish <= finish)
						{
							type.add(k->second.type);
//							k->second.valid = false;
							chrom.second.erase(k);
						}
					}
				}
			}


//			string & attName = i->second.tags[0].val;
			thisChromData.emplace(i->first,gtfRegion(i->second.start,finish,i->second.tags[0].val,i->second.strand,i->second.type));

		}
		
		//	

		//	And now for each region find the regions that it overlapsproduce overlap list:  The list of all regions that start before this region has ended.
		//	Then select the first of the regions, which will be used as the starting point for checking for region overlaps
		chromosomeEndIndexMap & thisChromEndMap = genomeEndIndex[chrom.first];

		//	We are using a temporary map of the address of the gtfEntries used to store the data.  This only works because the entries will not be moved
		//	during this process
		map<size_t,map<size_t,chromosomeGtfData::iterator> > tempMap;

		for (chromosomeGtfData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			//	Take the opportunity to produce a map of all the genes for holding counts

			if (!type.empty())
				geneCounts[i->second.name][type];

			//	And a parallel map of the ends of all of the gtfRegions/
			thisChromEndMap.emplace(i->second.finish,i);
		
			for (chromosomeGtfData::iterator j = next(i,1);(j != thisChromData.end()) && (j->first < i->second.finish);j++)
				tempMap[(size_t)&j->second].emplace(i->first,i);
		}
		for (chromosomeGtfData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			map<size_t,map<size_t,chromosomeGtfData::iterator> >::iterator j = tempMap.find((size_t)&i->second);

			if (j == tempMap.end())
				i->second.overlaps = i;
			else
				i->second.overlaps = j->second.begin()->second;
		}

	}

	geneCounts["__no_feature"][""];
	geneCounts["__ambiguous"][""];
	geneCounts["__too_low_aQual"][""];
	geneCounts["__not_aligned"][""];
	geneCounts["__alignment_not_unique"][""];

}

void gtfFileEx::outputChromData(const string & filename)
{
	TsvFile output;
	output.open(filename);

	for(auto i : genomeGtfData)
	{
		for (auto j : i.second)
			output.printEnd(i.first,j.first,j.second.finish,j.second.strand,j.second.name,j.second.type);
	}
}

