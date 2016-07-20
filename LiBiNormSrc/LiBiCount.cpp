#include <stdlib.h>
#include <crtdbg.h>
#include <ctime>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "LiBiCount.h"


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

bool gtfRegion::checkOverlap(const region & segment,bool & strict) const
{
	if ((segment.start >= start) && (segment.end <= finish))
	{
		strict = true;
		return true;
	}
	else if (((segment.start <= finish) && (segment.start >= start)) ||
				((segment.end <= finish) && (segment.end >= start)) ||
				((segment.start <= start) && (segment.end >= finish)))
	{
		strict = false;
		return true;
	}

	return false;
}

void regionLists::GetRegions(const BamAlignment & ba) {
	This[ba.RefID].GetRegions(ba);
}


void regionList::GetRegions(const BamAlignment & ba) {

	//	If we already have some regions then we need to combine them
	bool combine = size();

	//	Dont combine if the two reads are on reverse strands
	if (combine && (ba.IsReverseStrand() && (begin()->second.strand == '+')))
		combine = false;

	bool revStrand = ba.IsReverseStrand() != ba.IsFirstMate();

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
						combine = combineRegion(start,end-1,revStrand);
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
		combine = combineRegion(start,end-1,revStrand);
	else
		add(start,end-1,revStrand);
}



void gtfFileEx::index(mapZeroDef<string,size_t> & geneCounts)
{

	for (auto & chrom : entryMap)
	{
		//	In each chromosome
		chromosomeData & thisChromData = chromData[chrom.first];
		for (auto i = chrom.second.begin(); i != chrom.second.end();i++)
		{
			_DBG(size_t start = i->second.start;)
			size_t finish = i->second.finish;
			setEx<string> type(i->second.type);
			for (auto j = next(i,1);(j != chrom.second.end()) && (j->first < finish);)
			{
				auto k = j++;
				if (i->second.tags[0].val == k->second.tags[0].val)
				{
					if (k->second.finish > finish)
					{
						finish = k->second.finish;
						type.add(k->second.type);
						chrom.second.erase(k);
					}
					else if (k->second.finish <= finish)
					{
						type.add(k->second.type);
						chrom.second.erase(k);
					}
				}
			}


			string & attName = i->second.tags[0].val;
			thisChromData.emplace(i->first,gtfRegion(i->second.start,finish,i->second.tags[0].val,i->second.strand,move(type)));

		}
		
		//	

		//	And now for each region find the regions that it overlapsproduce overlap list:  The list of all regions that start before this region has ended.
		multimap<size_t,chromosomeData::iterator> & thisChromEndMap = chromEndIndex[chrom.first];

		for (chromosomeData::iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{

			//	Create entry if it does not exist
			geneCounts[i->second.name];

			thisChromEndMap.emplace(i->second.finish,i);
		
			for (chromosomeData::iterator j = next(i,1);(j != thisChromData.end()) && (j->first < i->second.finish);j++)
				j->second.overlaps.emplace(i->first,i);
		}
	}

	geneCounts["__no_feature"];
	geneCounts["__ambiguous"];
	geneCounts["__too_low_aQual"];
	geneCounts["__not_aligned"];
	geneCounts["__alignment_not_unique"];

}

void gtfFileEx::outputChromData(const string & filename)
{
	TsvFile output;
	output.open(filename);

	for(auto i : chromData)
	{
		for (auto j : i.second)
			output.printEnd(i.first,j.first,j.second.finish,j.second.strand,j.second.name,j.second.type);
	}
}

void LiBiCount::outputGeneCounts(const string & filename)
{
	TsvFile output;
	output.open(filename);

	for(auto i : geneCounts)
	{
		if (i.first.substr(0,2) != "__")
			output.printEnd(i.first,i.second);
	}
	geneCounts.print("__no_feature",output);
	geneCounts.print("__ambiguous",output);
	geneCounts.print("__too_low_aQual",output);
	geneCounts.print("__not_aligned",output);
	geneCounts.print("__alignment_not_unique",output);
}


void LiBiCount::addRead(const regionLists & segments,const gtfFileEx & gtfData)
{
	//	The paired end read consists of a number of segments.   If the two ends were aligned to different chromosomes
	//	then the segments will be on different chromosomes
	map<string,bool> genes;

	for (auto & chromSegments : segments)
	{
		//	For teh segments on each of the chromosomes (normally only one chromosome) get the gtfRegions for the chromosome
		map<string,chromosomeData>::const_iterator thisChromGtfRegions = gtfData.chromData.find(references[chromSegments.first].RefName);

		if (thisChromGtfRegions != gtfData.chromData.end())
		{

			//	Get the map of neds of gtfRegions associated with the chromosome
			const multimap<size_t,chromosomeData::iterator> & thisChromEndMap = gtfData.chromEndIndex.at(references[chromSegments.first].RefName);

			//	And find the first one that finishes at or beyond the start of the first segment
			multimap<size_t,chromosomeData::iterator>::const_iterator indirectIteratorStart = thisChromEndMap.lower_bound(chromSegments.second.begin()->second.start);

			//	And move back one to ensure we have the region that covers segment
			if (indirectIteratorStart != thisChromEndMap.begin())
				indirectIteratorStart--;

			chromosomeData::iterator gtfRegion = indirectIteratorStart->second;

			//	If this region overlaps any other regions then go to the one that starts the earliest.
			if (gtfRegion->second.overlaps.size())
				gtfRegion = gtfRegion->second.overlaps.begin()->second;


			//And now go through each of the segments
			for (auto & segment : chromSegments.second)
			{
				//	Trying out each of the gtfRegions in turn to see if there is an overlap
				auto j = gtfRegion;
				while ((j != thisChromGtfRegions->second.end()) && (j->first <= segment.second.end))
				{
					//	Check for strand match
					if (!useStrand || ((j->second.strand == segment.second.strand) != reverseStrand))
					{
						bool strict;
						if (j->second.checkOverlap(segment.second,strict))
						{
							//	If there are multiple matches and one is not strict then
							//	the overall match is not strict.
							auto g = genes.find(j->second.name);
							if (g != genes.end())
							{
								if (!strict)
									g->second = false;
							}
							else
								genes.emplace(j->second.name,strict);
							gtfRegion = j;
						}
					}
					j++;
				}
			}
		}
	}

	switch(countMode)
	{
	case intersect_union:
		if (genes.size() == 1)
		{
			_DBG(
				if (genes.begin()->first == "TLE4")
				cout << regions.name << endl;)
			geneCounts[genes.begin()->first]++;
		}
		else if (genes.size() > 1)
			geneCounts["__ambiguous"]++;
		else
			geneCounts["__no_feature"]++;
		break;
	case intersect_strict:
	case intersect_nonempty:
		{
			vector<string> strictNames,nonemptyNames;
			for (auto & gene : genes)
			{
				if (gene.second)
					strictNames.push_back(gene.first);
				else if (gene.second)
					nonemptyNames.push_back(gene.first);
			}
			if (strictNames.size() == 1)
				geneCounts[strictNames[0]]++;
			else if ((nonemptyNames.size() ==1) && (countMode == intersect_nonempty))
				geneCounts[nonemptyNames[0]]++;
			break;
		}
	}
}

int LiBiCount::main(int argc, char **argv)
{
	stringEx bamFileName,gtfFileName,
		id_attribute = "gene_name";

	setEx<string> feature_type("exon");

	reverseStrand = false;
	useStrand = true;
	countMode = intersect_union;
//	setEx<string> feature_type("exon","five_prime_utr","three_prime_utr");

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if(argc == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm:    RNA-seq library bias normalisation   \n");
		return EXIT_SUCCESS;
	}

	int ni = 1;
	while(ni < argc)
	{
		if(strcmp(argv[ni], "-b") == 0)
		{
			bamFileName = argv[++ni];
		}
		else if(strcmp(argv[ni], "-g") == 0)
		{
			gtfFileName = argv[++ni];
		}
		else
		{
			cout << "Invalid parameter";
		}
		ni++;
	}

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();

	clock_t begin = clock();

	if (!genomeDef.open(gtfFileName,id_attribute,feature_type))
		exitFail("Could not open gtf file: ",gtfFileName);

	genomeDef.index(geneCounts);

	cout << "GFF file sonsolidated." << endl;

	clock_t now = clock();
	double elapsed_secs = double(now - begin) / CLOCKS_PER_SEC;
	cout << "Elapsed time " << elapsed_secs << endl;

//	genomeDef.outputChromData(gtfFileName.replaceSuffix(".txt"));


	processBamData();


	outputGeneCounts(bamFileName.replaceSuffix(".counts.txt"));

	now = clock();
	elapsed_secs = double(now - begin) / CLOCKS_PER_SEC;

	cout << "Elapsed time " << elapsed_secs << endl;

	string test;
	cin >> test;

	return EXIT_SUCCESS;
}


void LiBiCount::processBamData()
{
	BamAlignment ba1;
	BamAlignment ba2;

	size_t samCounter(0);

	bool OK = reader.GetNextAlignment(ba1);

	while (OK)
	{
		_DBG(string name = ba1.Name;
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1102:10307:46089");)

			
		regionLists regions;
		_DBG(regions.name = name;)

		bool readAlreadyRead = false;
		
		regions.GetRegions(ba1);

		if (ba1.IsFirstMate())
		{
			reader.GetNextAlignment(ba2);

			if (ba2.Name == ba1.Name)
			{
				regions.GetRegions(ba2);
			}
			else
			{
				readAlreadyRead = true;
			}
		}

		if ((++samCounter % 100000) == 0)
			cout << samCounter << " SAM alignment record pairs processed." << endl;

		addRead(regions,genomeDef);

		if (readAlreadyRead)
			ba1 = ba2;
		else
			OK = reader.GetNextAlignment(ba1);
	}

}

