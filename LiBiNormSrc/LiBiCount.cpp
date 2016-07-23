#include <stdlib.h>
#include <crtdbg.h>
#include "libCommon.h"
#include "stringEx.h"
#include "containerEx.h"
#include "LiBiCount.h"


int LiBiCount::main(int argc, char **argv)
{
	string test;

	stringEx bamFileName,gtfFileName,outputFilename,
		id_attribute = "gene_name";

	setEx<string> feature_type("exon");

	reverseStrand = true;
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
		if(strcmp(argv[ni], "-f") == 0)
		{
			fileCompare();
			exitSuccess();
		}
		else if(strcmp(argv[ni], "-b") == 0)
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
		else if((strcmp(argv[ni], "-o") == 0) || (strcmp(argv[ni], "--samout") == 0))
		{
			outputFilename = argv[++ni];
		}

		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}

	if (outputFilename && !outputFile.open(outputFilename))
		exitFail("Unable to open output file: ",outputFilename);

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();

	initClock();

//	clock_t begin = clock();

	if (!genomeDef.open(gtfFileName,verbose,id_attribute,feature_type))
		exitFail("Could not open gtf file: ",gtfFileName);

	genomeDef.index(geneCounts);

//	_DBG(cin >> test;)

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
		size_t length;
		overlapCounts(size_t partial=0,size_t strict=0,size_t length = 0):partial(partial),strict(strict),length(length){};
	};

	struct chromosomeGeneInfo: public map<string,overlapCounts>
	{
		size_t noMatch;
		size_t nSegments;
		chromosomeGeneInfo():noMatch(0),nSegments(0){};
	} genes[2];

	bool standardPair = true;
	bool nonOverlappingGenes = false;

	size_t pairNo = 0;
	for (auto & chromSegments : segments)
	{
		if (!chromSegments.second.standardPair)
			standardPair = false;
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
				genes[pairNo].nSegments++;
				vector<gtfOverlap> overlaps;
				//	Trying out each of the gtfRegions in turn to see if there is an overlap
				auto j = gtfRegion;
				while ((j != thisChromGtfRegions->second.end()) && (j->first <= segment.second.end))
				{
					//	Check for strand match
					if (!useStrand || ((j->second.strand == segment.second.strand) != reverseStrand))
					{
						bool strict;
						if (j->second.checkOverlap(segment.second,strict,overlaps))
						{
							//	If there are multiple matches and one is not strict then
							//	the overall match is not strict.
/*							auto g = genes[pairNo].find(j->second.name);
							if (g != genes[pairNo].end())
							{
								g->second.partial++;
								if (strict)
									g->second.strict++;
							}
							else
								genes[pairNo].emplace(j->second.name,overlapCounts(1,strict?1:0));*/
							if (j->second.start > segment.second.end)
								gtfRegion = j;
						}
					}
					j++;
				}
				if (overlaps.size() == 0)
				{
					genes[pairNo].noMatch++;
				}
/*				else if (overlaps.size() == 1)
				{
					size_t length = overlaps[0].gtfReg.finish - overlaps[0].gtfReg.start;
					overlapCounts & olc = genes[pairNo][overlaps[0].gtfReg.name];
					olc.partial++;
					olc.strict += overlaps[0].strict?1:0;
					olc.length += length;*/

/*					auto g = genes[pairNo].find(overlaps[0].gtfReg.name);
					size_t length = overlaps[0].gtfReg.finish - overlaps[0].gtfReg.start;
					if (g == genes[pairNo].end())
					{
						genes[pairNo].emplace(overlaps[0].gtfReg.name,overlapCounts(1,overlaps[0].strict?1:0,length));
					}
					else
					{
						g->second.partial++;
						g->second.strict += overlaps[0].strict?1:0;
						g->second.length += length;
					}*/
/*				}
				else*/
				{
					size_t min = INT_MAX,max = 0;
					stringEx gene;
	
					for (size_t x = 0;x < overlaps.size();x++)
					{
						//	Create dummy entry for every gene;
						genes[pairNo][overlaps[x].gtfReg.name];

						if (overlaps[x].strict)
							genes[pairNo][overlaps[x].gtfReg.name].strict++;


						if ((overlaps[x].start == min) && (overlaps[x].finish == max))
						{
							//Identical
							gene = "";

						}
						else if ((overlaps[x].start >= min) && (overlaps[x].finish <= max))
						{
							//	smaller ignore
						}
						else if ((overlaps[x].start <= min) && (overlaps[x].finish >= max))
						{
							//	bigger replace
							min = overlaps[x].start,max = overlaps[x].finish;
							gene = overlaps[x].gtfReg.name;
						}
						else
						{
							nonOverlappingGenes = true;
						}
					}
					if (gene)
					{
						overlapCounts & olc = genes[pairNo][gene];
						olc.partial++;
						olc.length += (max-min);
/*
						auto g = genes[pairNo].find(gene);
						if (g == genes[pairNo].end())
						{
							genes[pairNo].emplace(gene,overlapCounts(1,strict?1:0,max-min));
						}
						else
						{
							g->second.partial++;
							g->second.strict += strict?1:0;
							g->second.length += (max-min);
						}*/
//						if (g != genes[pairNo].end())
//							g->second.length = max - min;
					}
				}
			}
		}
		pairNo++;
	}


	switch(countMode)
	{
	case intersect_strict:
	case intersect_nonempty:		
		{
			if ((countMode == intersect_nonempty) || standardPair)
			{
				mapZeroDef <string,size_t> strictNames;

				size_t totSegments = genes[0].nSegments + genes[1].nSegments;

				for (size_t chromosome = 0;chromosome < 2;chromosome++)
				{
					for (auto & gene : genes[chromosome])
					{
						if (gene.second.strict == genes[chromosome].nSegments)
							strictNames[gene.first] += gene.second.strict;
					}
				}
				for (auto i = strictNames.begin(); i != strictNames.end();)
				{
					if (i->second != totSegments)
					{
						auto j = i++;
						strictNames.erase(j);
					}
					else
						i++;
				}

				if (strictNames.size() == 1)
				{
					const string & name = strictNames.begin()->first;
					geneCounts[name]++;
					if (outputFile.is_open())
						outputFile.printEnd("strict",name,segments.name);
					break;
				}
				else if (strictNames.size() > 1)
				{
					if (outputFile.is_open())
						outputFile.printEnd("strict","__ambiguous",segments.name);
					geneCounts["__ambiguous"]++;
					break;
				}
				else if (countMode == intersect_strict)
				{
					if (outputFile.is_open())
						outputFile.printEnd("strict","__no_feature",segments.name);
					geneCounts["__no_feature"]++;	
					break;
				}
			}
			//  If we have not found a match, and it is not strict, run on and try union
		}
	case intersect_union:
		{
			size_t Ngenes = genes[0].size() + genes[1].size();
			if ((Ngenes == 1) && (standardPair || (countMode == intersect_nonempty)))
			{
				size_t readNo = genes[1].size();
				if (outputFile.is_open())
					outputFile.printEnd("strict",genes[readNo].begin()->first,segments.name);
				geneCounts[genes[readNo].begin()->first]++;
			}
			else if (Ngenes > 1)
			{
				if (countMode == intersect_union)
				{
					if (outputFile.is_open())
						outputFile.printEnd("union","__ambiguous",segments.name);
					geneCounts["__ambiguous"]++;
				}
				else 
				{
					size_t bestLength = 0;
					stringEx bestGene;
					int maxLength = 0;
					bool ambiguous = false;

					if ((genes[0].size() == 0) || (genes[1].size() == 0))
					{
						size_t chromId = (genes[0].size())?0:1;

						for (auto & gene : genes[chromId])
						{
							if (gene.second.length > maxLength)
							{
								if (maxLength > 0)	
									ambiguous = true;
								else
									maxLength = gene.second.length;
							}
							if (gene.second.partial == (genes[chromId].nSegments - genes[chromId].noMatch))
							{
								if (gene.second.length > bestLength)
								{
									bestGene = gene.first;
									bestLength = gene.second.length;
								}
								else if (gene.second.length == bestLength)
								{
									ambiguous = true;
								}
							}
						}
					}
					if (ambiguous || nonOverlappingGenes)
					{
						if (outputFile.is_open())
							outputFile.printEnd("non_empty","__ambiguous",segments.name);
						geneCounts["__ambiguous"]++;
					}
					else if (bestGene)
					{
						if (outputFile.is_open())
							outputFile.printEnd("non_empty",bestGene,segments.name);
						geneCounts[bestGene]++;
					}
					else
					{
						if (outputFile.is_open())
							outputFile.printEnd("non_empty","__no_feature",segments.name);
						geneCounts["__no_feature"]++;
					}
				}
			}
			else
			{
				if (outputFile.is_open())
					outputFile.printEnd("union","__no_feature",segments.name);
				geneCounts["__no_feature"]++;
			}
			break;
		}
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
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1101:9276:56896");)

			
		regionLists regions;
		regions.name = ba1.Name;

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



void LiBiCount::fileCompare()
{
	ifstream samFile;
	samFile.open("Y:\\SysmedIBD\\CD\\test\\nonemptyTest\\nonempty_large.sam");
	if (!samFile.is_open()) return;

	ifstream myFile;
	myFile.open("Y:\\SysmedIBD\\CD\\test\\nonemptyTest\\matches_2.txt");
	if (!myFile.is_open()) return;


	TsvFile outFile;
	outFile.open("Y:\\SysmedIBD\\CD\\test\\nonemptyTest\\comparison_2.txt");

	string line;

	vectorEx<size_t> cols(1,2,3,4,5);


	vector<string> myParams;
	vector<stringEx> samParams[2];

	bool OK = getline(samFile,line);
	parseTsv(line,samParams[0]);

	while (OK) {

		getline(myFile,line);
		myParams.clear();
		parseTsv(line,myParams);


		getline(samFile,line);
		samParams[1].clear();
		parseTsv(line,samParams[1]);

		if (samParams[0][0] != myParams[2])
			exitFail("Mismatch ",samParams[0][0],myParams[2]);

		outFile.printMiddle(myParams[2],myParams[0],myParams[1]);
		for (size_t i : cols)
		{
				outFile.printMiddle(samParams[0][i]);
		}
		for (size_t i = 8;i < samParams[0].size();i++)
		{
			if (samParams[0][i].startsWith("XF:Z:"))
			{
				outFile.printMiddle(samParams[0][i].substr(5));
			}
		}
		outFile.printEnd();

		if (samParams[1][0] != myParams[2])
		{
			swap(samParams[0],samParams[1]);
		}
		else
		{
			OK = getline(samFile,line);
			samParams[0].clear();
			parseTsv(line,samParams[0]);
		}

	};

}
