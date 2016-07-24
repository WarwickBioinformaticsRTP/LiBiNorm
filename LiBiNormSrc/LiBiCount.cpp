#include <stdlib.h>

#ifdef _WIN32
#include <crtdbg.h>
#endif

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
			fileCompare(argv[++ni]);
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
	} genes;

	bool nonOverlappingGenes = false;

	size_t pairNo = 0;
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

				genes.nSegments++;
				vector<gtfOverlap> overlaps;
				//	Trying out each of the gtfRegions in turn to see if there is an overlap

				for (auto j = gtfRegion; (j != thisChromGtfRegions->second.end()) && (j->first <= segment.second.end); j++)
				{
					//	Check for strand match
					if (!useStrand || ((j->second.strand == segment.second.strand) != reverseStrand))
					{
						j->second.checkOverlap(segment.second,overlaps);
					}

					gtfRegion = j->second.overlaps;
				}

				if (overlaps.size() == 0)
				{
					genes.noMatch++;
				}
				else
				{
					size_t min = INT_MAX,max = 0;
					set<string> geneSet;
	
					for (auto & overlap: overlaps)
					{
						//	Create dummy entry for every gene;
						genes[overlap.geneName];

						if (overlap.strict)
							genes[overlap.geneName].strict++;


						if ((overlap.start == min) && (overlap.finish == max))
						{
							//Identical
							geneSet.emplace(overlap.geneName);
						}
						else if ((overlap.start >= min) && (overlap.finish <= max))
						{
							//	smaller ignore
						}
						else if ((overlap.start <= min) && (overlap.finish >= max))
						{
							//	bigger replace
							min = overlap.start;
							max = overlap.finish;
							geneSet.clear();
							geneSet.emplace(overlap.geneName);
						}
						else
						{
							nonOverlappingGenes = true;
						}
					}
					if (geneSet.size())
					{
						for (auto gene:geneSet)
						{
							genes[gene].partial++;
							genes[gene].length += (max-min);
						}
					}
				}
			}
		}
		pairNo++;
	}

	static string ambiguousString = "__ambiguous";
	static string noFeatureString = "__no_feature";

	switch(countMode)
	{
	case intersect_strict:
	case intersect_nonempty:		
		{

			const string * result = nullptr;

			for (auto & gene : genes)
			{
				if (gene.second.strict == genes.nSegments)
				{
					if (result)
					{
						result = &ambiguousString;
						break;
					}
					else
						result = &gene.first;
				}
			}


			if ((countMode == intersect_strict) && (result == nullptr))
				result = &noFeatureString;

			if (result != nullptr)
			{
				geneCounts[*result]++;
				if (outputFile.is_open())
					outputFile.printEnd("strict",*result,segments.name);
				break;
			}

		}
	case intersect_union:
		{
			size_t Ngenes = genes.size();
			const string * result = nullptr;

			if (Ngenes == 1)
			{
				result = &genes.begin()->first;
			}
			else if (Ngenes > 1)
			{
				if (countMode == intersect_union)
				{
					result = &ambiguousString;
				}
				else 
				{
					size_t bestLength = 0;

					result = &noFeatureString;

					for (auto & gene : genes)
					{
						if (gene.second.partial == (genes.nSegments - genes.noMatch))
						{
							if (gene.second.length > bestLength)
							{
								result = &gene.first;
								bestLength = gene.second.length;
							}
							else if (gene.second.length == bestLength)
							{
								result = &ambiguousString;
							}
						}
					}

					if (nonOverlappingGenes && (result != &ambiguousString))
					{
						result = &noFeatureString;
					}


					if (outputFile.is_open())
						outputFile.printEnd("non_empty",*result,segments.name);
					geneCounts[*result]++;
					break;
				}
			}
			else
			{
				result = &noFeatureString;
			}

			if (outputFile.is_open())
				outputFile.printEnd("union",*result,segments.name);
			geneCounts[*result]++;
			break;
		}
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
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1102:5557:44071");)

			
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
			cout << samCounter << " BAM alignment record pairs processed." << endl;

		addRead(regions,genomeDef);

		if (readAlreadyRead)
			swap(ba1,ba2);
		else
			OK = reader.GetNextAlignment(ba1);
	}

}



void LiBiCount::fileCompare(const string & mode)
{
	ifstream samFile;
	stringEx filename("Y:\\SysmedIBD\\CD\\test\\",mode,"Test\\",mode,"_large.sam");
	samFile.open(filename);
	if (!samFile.is_open())
		exitFail("Failed to open samfile ",filename);

	ifstream myFile;
	myFile.open(stringEx("Y:\\SysmedIBD\\CD\\test\\",mode,"Test\\matches_2.txt"));
	if (!myFile.is_open())
		exitFail("Failed to txt samfile");



	TsvFile outFile;
	outFile.open(stringEx("Y:\\SysmedIBD\\CD\\test\\",mode,"Test\\comparison_2.txt"));

	string line;

	vectorEx<size_t> cols(1,2,3,4,5);


	vector<string> myParams;
	vector<stringEx> samParams[2];

	bool OK = getline(samFile,line).good();

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
			OK = getline(samFile,line).good();
			samParams[0].clear();
			parseTsv(line,samParams[0]);
		}

	};

}
