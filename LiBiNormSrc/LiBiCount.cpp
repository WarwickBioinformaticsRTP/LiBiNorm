#include <stdlib.h>
#include <crtdbg.h>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "LiBiCount.h"


int LiBiCount::main(int argc, char **argv)
{
	stringEx bamFileName,gtfFileName,
		id_attribute = "gene_name";

	setEx<string> feature_type("exon");

	reverseStrand = false;
	useStrand = true;
	verbose = true;
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
		else if((strcmp(argv[ni], "-m") == 0) || (strcmp(argv[ni], "--mode") == 0))
		{
			string mode = argv[++ni];
			if (mode == "union")
				countMode = intersect_union;
			else if (mode == "intersection-strict")
				countMode = intersect_strict;
			else if (mode == "intersection-nonempty")
				countMode = intersect_nonempty;
			else
				exitFail("Invalid mode: ",mode);
		}
		else if((strcmp(argv[ni], "-q") == 0) || (strcmp(argv[ni], "--quiet") == 0))
		{
			verbose = false;
		}
		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();

	initClock();

//	clock_t begin = clock();

	if (!genomeDef.open(gtfFileName,verbose,id_attribute,feature_type))
		exitFail("Could not open gtf file: ",gtfFileName);

	genomeDef.index(geneCounts);

	if (verbose)
	{
		cerr << "GFF file sonsolidated." << endl;
		elapsedTime();
	}

	_DBG(genomeDef.outputChromData(gtfFileName.replaceSuffix(".txt"));)

	processBamData();

	outputGeneCounts(bamFileName.replaceSuffix(".counts.txt"));

	if (verbose)
		elapsedTime();

	string test;
	_DBG(cin >> test;)

	return EXIT_SUCCESS;
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

	struct overlapCounts
	{
		size_t partial;
		size_t strict;
		overlapCounts(size_t partial,size_t strict):partial(partial),strict(strict){};
	};

	map<string,overlapCounts> genes;
	size_t Nsegments = 0;

	for (auto & chromSegments : segments)
	{
		//	For teh segments on each of the chromosomes (normally only one chromosome) get the gtfRegions for the chromosome
		genomeGtfRegions::const_iterator thisChromGtfRegions = gtfData.genomeGtfData.find(references[chromSegments.first].RefName);

		if (thisChromGtfRegions != gtfData.genomeGtfData.end())
		{

			//	Get the map of neds of gtfRegions associated with the chromosome
			const chromosomeEndIndexMap & thisChromEndMap = gtfData.genomeEndIndex.at(references[chromSegments.first].RefName);

			//	And find the first one that finishes at or beyond the start of the first segment
			chromosomeEndIndexMap::const_iterator indirectIteratorStart = thisChromEndMap.lower_bound(chromSegments.second.begin()->second.start);

			//	And move back one to ensure we have the region that covers segment
			if (indirectIteratorStart != thisChromEndMap.begin())
				indirectIteratorStart--;

			chromosomeGtfData::iterator gtfRegion = indirectIteratorStart->second;

			//	If this region overlaps any other regions then go to the one that starts the earliest.
			//	If there were no overlaps then default is the overlaps points to self
			gtfRegion = gtfRegion->second.overlaps;


			//And now go through each of the segments
			for (auto & segment : chromSegments.second)
			{
				Nsegments++;
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
								g->second.partial++;
								if (strict)
									g->second.strict++;
							}
							else
								genes.emplace(j->second.name,overlapCounts(1,strict?1:0));
							if (j->second.start > segment.second.end)
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
	case intersect_strict:
	case intersect_nonempty:		
		{
			vector<string> strictNames;
			for (auto & gene : genes)
			{
				if (gene.second.strict == Nsegments)
					strictNames.push_back(gene.first);
			}
			if (strictNames.size() == 1)
			{
				geneCounts[strictNames[0]]++;
				break;
			}
			else if (strictNames.size() > 1)
			{
				geneCounts["__ambiguous"]++;
				break;
			}
			else if (countMode == intersect_strict)
			{
				geneCounts["__no_feature"]++;	
				break;
			}
			//  If we have not found a match, and it is not strict, run on and try union
		}
	case intersect_union:
		if (genes.size() == 1)
			geneCounts[genes.begin()->first]++;
		else if (genes.size() > 1)
		{
			if (countMode == intersect_union)
				geneCounts["__ambiguous"]++;
			else 
			{
/*				set<string> partialNames;
				for (auto & gene : genes)
				{
					if (gene.second.partial == Nsegments)
						partialNames.emplace(gene.first);
				}
				if (partialNames.size() > 1)
					geneCounts["__ambiguous"]++;
				else if (partialNames.size() == 1)
					geneCounts[*partialNames.begin()]++;
				else*/
					geneCounts["__ambiguous"]++;
			}
		}
		else
			geneCounts["__no_feature"]++;
		break;
/*	case intersect_nonempty:
		{
			vector<string> strictNames,nonemptyNames;
			for (auto & gene : genes)
			{
				if (gene.second == Nsegments)
					strictNames.push_back(gene.first);
				else
					nonemptyNames.push_back(gene.first);
			}
			if (strictNames.size() == 1)
				geneCounts[strictNames[0]]++;
			else if (nonemptyNames.size() ==1)
				geneCounts[nonemptyNames[0]]++;
			else if ((strictNames.size() > 1) || ((strictNames.size() == 0 ) && (nonemptyNames.size() > 1)))
				geneCounts["__ambiguous"]++;
			else
				geneCounts["__no_feature"]++;
			break;
		}*/
	}
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
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1101:18840:53081");)

			
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
			swap(ba1,ba2);
		else
			OK = reader.GetNextAlignment(ba1);
	}

}

