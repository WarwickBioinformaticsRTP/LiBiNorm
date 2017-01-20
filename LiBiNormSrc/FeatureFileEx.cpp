#include "FeatureFileEx.h"


using namespace std;

void featureRegion::checkOverlap(const region & segment,vector<featureOverlap> & overlapList) const
{
	rna_pos_type RNAstartPos = max<long>(segment.start - start + RNAstart,0);
	rna_pos_type RNAendPos = segment.end-start + RNAstart;

	//	strict fit
	if ((segment.start >= start) && (segment.end <= finish))
	{
		overlapList.emplace_back(segment.start,segment.end, RNAstartPos, RNAendPos,true,name,type);
	}
	//	non-strict fits
	else if ((segment.start <= finish) && (segment.start >= start)) 
	{
		overlapList.emplace_back(segment.start,finish, RNAstartPos, RNAendPos, false,name,type);
	}
	else if ((segment.end <= finish) && (segment.end >= start))
	{
		overlapList.emplace_back(start,segment.end, RNAstartPos, RNAendPos, false,name,type);
	}
	else if ((segment.start <= start) && (segment.end >= finish))
	{
		overlapList.emplace_back(start,finish, RNAstartPos, RNAendPos,false,name,type);
	}
}

featureRegion::featureRegion(size_t start, size_t finish, const std::string & name, char strand, const std::string & type) :
	start(start), finish(finish),  RNAstart(0), name(name), type(type), strand(strand)
{
	overlaps = new chromosomeFeatureData::iterator();
};
featureRegion::~featureRegion() {
	delete (overlaps);
};


void featureFileEx::index(GeneCountData & geneCounts)
{
	//	If there are already entries in the geneCounts data at this stage it is because
	//	we have preloaded them with a set of genes/transcripts that we are specifically
	//	interested in.  At this point we then get rid of the rest
	if (geneCounts.readPositionData.size() > 1)
	{
		for (auto & chrom : entryMap)
		{
			{
				//	First get rid of entries associated with genes that we are not interested in 
				for (auto i = chrom.second.begin(); i != chrom.second.end();)
				{
					if (geneCounts.readPositionData.find(i->second.tags[0].val) == geneCounts.readPositionData.end())
					{
						//			if (!geneList.contains(i->second.tags[0].val))
						auto j = i++;
						chrom.second.erase(j);
					}
					else
						i++;
				}
			}
		}
	}



	for (auto & chrom : entryMap)
	{
		//	In each chromosome go through all of the regions to see what regions can be amalgamated
		chromosomeFeatureData & thisChromData = genomeGtfData[chrom.first];
		for (auto i = chrom.second.begin(); i != chrom.second.end();i++)
		{

			size_t finish = i->second.finish;
			setEx<string> type{{i->second.type}};
			for (auto j = next(i,1);(j != chrom.second.end()) && (j->first <= (finish + 1));)
			{
				auto k = j++;
				if (k != chrom.second.end())
				{
					//	It turns out that the Yam1 gene is defined on both strands, so need to check strands
					if ((i->second.type == k->second.type) && (i->second.tags[0].val == k->second.tags[0].val) && (i->second.strand == k->second.strand))
					{
						//	If the second region extends beyond the firat then make the region longer and add
						//	the second type to the list of types associated with the region
						if (k->second.finish > finish)
						{
							finish = k->second.finish;
							type.add(k->second.type);
							chrom.second.erase(k);
						}
						//	If teh second region is the same length or shorter then just add the type associated with the new region
						else if (k->second.finish <= finish)
						{
							type.add(k->second.type);
							chrom.second.erase(k);
						}
					}
				}
			}
			thisChromData.emplace(i->first,featureRegion(i->second.start,finish,i->second.tags[0].val,i->second.strand,i->second.type));

		}
		
		//	

		//	And now for each region find the regions that it overlaps and produce overlap list:  The list of all regions that start before this region has ended.
		//	Then select the first of the regions, which will be used as the starting point for checking for region overlaps
		chromosomeEndIndexMap & thisChromEndMap = genomeEndIndex[chrom.first];

		
		
		//	We are using a temporary map of the address of the gtfEntries used to store the data.  This only works because the entries will not be moved
		//	during this process
		map<void *,map<size_t,chromosomeFeatureData::iterator> > tempMap;

		for (chromosomeFeatureData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{

			//	Take the opportunity to produce a map of all the genes for holding counts
			geneCounts.addEntry(i->second.name);

			//	And a parallel map of the ends of all of the featureRegions/
			thisChromEndMap.emplace(i->second.finish,i);

		
			for (chromosomeFeatureData::iterator j = next(i,1);(j != thisChromData.end()) && (j->first < i->second.finish);j++)
				tempMap[&j->second].emplace(i->first,i);
		}
		for (chromosomeFeatureData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			VEC_DATA_TYPE & length = geneCounts.lengths[0].at(geneCounts.readPositionData.at(i->second.name).index);
			auto j = tempMap.find(&i->second);

			if (j == tempMap.end())
			{
				*i->second.overlaps = i;
				//	Add featureRegion to the list of regions associated with the gene
				genes[i->second.name].addRegion(&i->second, false,length);
			}
			else
			{
				*i->second.overlaps = j->second.begin()->second;
				genes[i->second.name].addRegion(&i->second, true,length);
			}
		}

	}
	geneCounts.addErrorEntry(noFeatureString);
	geneCounts.addErrorEntry(ambiguousString);
	geneCounts.addErrorEntry(lowQualString);
	geneCounts.addErrorEntry(notAlignedString);
	geneCounts.addErrorEntry(notUnique);

}

void featureFileEx::outputChromData(const string & filename)
{
	TsvFile output;
	output.open(filename);

	if (!output.is_open())
		exitFail("Unable to open output file",filename);

	for(auto & i : genomeGtfData)
	{
		for (auto & j : i.second)
			output.printEnd(i.first,j.first,j.second.finish,j.second.strand,j.second.name,j.second.type);
	}
}

