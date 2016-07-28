#include <stdlib.h>

#ifdef _WIN32
#include <crtdbg.h>
#endif

#include "libCommon.h"
#include "LiBiCount.h"
#include "parser.h"

#ifdef _DEBUG
#define READ_CACHE_SIZE 100000
#else
#define READ_CACHE_SIZE 2000000
#endif

using namespace std;

int LiBiCount::main(int argc, char **argv)
{
	std::string test;

	stringEx bamFileName,gtfFileName,outputFilename,
		id_attribute = "gene_name";

	setEx<string> feature_type;

	reverseStrand = true;
	useStrand = true;
	verbose = true;
	countMode = intersect_union;
	maxCacheSize = READ_CACHE_SIZE;

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
		if(strcmp(argv[ni], "-t") == 0)
		{
			fileCompare(argv[++ni]);
			exitSuccess();
		}
		else if(strcmp(argv[ni], "-c") == 0)
		{
			maxCacheSize = atoi(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-f") == 0)
		{
			feature_type.emplace(argv[++ni]);
		}
		else if(strcmp(argv[ni], "-i") == 0)
		{
			id_attribute = argv[++ni];
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
		else if((strcmp(argv[ni], "-r") == 0) || (strcmp(argv[ni], "--results") == 0))
		{
			resultsFilename = argv[++ni];
		}

		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	if (feature_type.size() == 0)
		feature_type.emplace("exon");

	if (!resultsFilename)
		exitFail("No results filename specified.");

	if (outputFilename && !outputFile.open(outputFilename))
		exitFail("Unable to open output file: ",outputFilename);

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();

	initClock();

	if (!genomeDef.open(gtfFileName,verbose,id_attribute,feature_type))
		exitFail("Could not open gtf file: ",gtfFileName);

	genomeDef.index(geneCounts,(feature_type.size()==1)?*feature_type.begin():"");

//	_DBG(cin >> test;)

	if (verbose)
	{
		cerr << "GFF file sonsolidated." << endl;
		elapsedTime();
	}

	_DBG(genomeDef.outputChromData(gtfFileName.replaceSuffix(".txt"));)

	if (!processOrderedBamData())
	{
		reader.Rewind();
		geneCounts.reset();
		if (verbose)
			cerr << "Processing data assuming that it is not ordered by read name." << endl;
		processUnorderedBamData();
	}

	if(!outputGeneCounts(resultsFilename))
		exitFail("Unable to output counts to :",resultsFilename);

	if (verbose)
		elapsedTime();

	cin >> test;

	return EXIT_SUCCESS;
}



bool LiBiCount::outputGeneCounts(const string & filename)
{
	TsvFile output;
	
	if (!output.open(filename))
		return false;

	for(auto i : geneCounts)
	{
		for(auto j : i.second)
		{
			if (i.first.substr(0,2) != "__")
				output.printEnd(i.first,j.first,j.second);
		}
	}
	geneCounts.print("__no_feature",output);
	geneCounts.print("__ambiguous",output);
	geneCounts.print("__too_low_aQual",output);
	geneCounts.print("__not_aligned",output);
	geneCounts.print("__alignment_not_unique",output);
	return true;
}


//	overlapCounts and chromosomeInfo would have been declared inside addRead as they are local to addRead.  However this produces a C++ warning that the
//	decorated name is too long, which can cause problems with debugging

struct overlapCounts
{
	size_t partial;
	size_t strict;
	size_t length;
	overlapCounts(size_t partial=0,size_t strict=0,size_t length = 0):partial(partial),strict(strict),length(length){};
};

struct chromosomeGeneInfo: public map<string,map <string,overlapCounts> >
{
	size_t noMatch;
	size_t nSegments;
	chromosomeGeneInfo():noMatch(0),nSegments(0){};
	size_t size()
	{
		size_t _s = 0;
		for (auto & i : This)
		{
			_s += i.second.size();
		}
		return _s;
	}
};


void LiBiCount::addRead(const regionLists & segments,const gtfFileEx & gtfData)
{
	//	The paired end read consists of a number of segments.   If the two ends were aligned to different chromosomes
	//	then the segments will be on different chromosomes
	chromosomeGeneInfo genes;

	bool nonOverlappingGenes = false;

	size_t pairNo = 0;
	for (auto & chromSegments : segments.data)
	{
		//	For teh segments on each of the chromosomes (normally only one chromosome) get the gtfRegions for the chromosome
		genomeGtfRegions::const_iterator thisChromGtfRegions = gtfData.genomeGtfData.find(references[chromSegments.first].RefName);

		if (thisChromGtfRegions != gtfData.genomeGtfData.end())
		{

			//	Get the map of neds of gtfRegions associated with the chromosome
			const chromosomeEndIndexMap & thisChromEndMap = gtfData.genomeEndIndex.at(references[chromSegments.first].RefName);

			//	And find the first one that finishes at or beyond the start of the first segment
			chromosomeEndIndexMap::const_iterator indirectIteratorStart = thisChromEndMap.lower_bound(chromSegments.second.data.begin()->second.start);

			//	And move back one to ensure we have the region that covers segment
			if (indirectIteratorStart != thisChromEndMap.begin())
				indirectIteratorStart--;

			chromosomeGtfData::iterator gtfRegion = indirectIteratorStart->second;

			//	If this region overlaps any other regions then go to the one that starts the earliest.
			//	If there were no overlaps then default is the overlaps points to self
			gtfRegion = gtfRegion->second.overlaps;


			//And now go through each of the segments
			for (auto & segment : chromSegments.second.data)
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
					set<pair<string,string> > geneSet;
	
					for (auto & overlap: overlaps)
					{
						//	Create dummy entry for every region type that the read overlapped;
						genes[overlap.geneName][overlap.type];

						if (overlap.strict)
							genes[overlap.geneName][overlap.type].strict++;


						if ((overlap.start == min) && (overlap.finish == max))
						{
							//Identical
							geneSet.emplace(overlap.geneName,overlap.type);
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
							geneSet.emplace(overlap.geneName,overlap.type);
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
							genes[gene.first][gene.second].partial++;
							genes[gene.first][gene.second].length += (max-min);
						}
					}
				}
			}
		}
		pairNo++;
	}

	static string blankString = "";

	//	Result options
	static string ambiguousString = "__ambiguous";
	static string noFeatureString = "__no_feature";

	//	Result mode options
	static string strictString = "strict";
	static string nonemptyString = "nonempty";
	static string unionString = "union";

	const string * result = nullptr;
	const string * type = &blankString;
	const string * mode = &blankString;

	switch(countMode)
	{
	case intersect_strict:
	case intersect_nonempty:		
		{

			//	For a strict match, all of the read segments must lie inside an annotated region of the same gene
			for (auto & gene : genes)
			{
				for (auto & regionType : gene.second)
				{
					if (regionType.second.strict == genes.nSegments)
					{
						if (result)
						{
							//	If we have two strict matches then the result is ambigous, no need to look any further
							result = &ambiguousString;
							type = &blankString;
							break;
						}
						else
						{
							//	A strict match, keep looking as there may be more
							result = &gene.first;
							type = &regionType.first;
						}
					}
				}
			}

			//  If this is strict mode then the only option is a strict match
			if ((countMode == intersect_strict) && (result == nullptr))
				result = &noFeatureString;

			if (result != nullptr)
			{
				mode = &strictString;
				break;
			}

			//	Note that if the mode is intersection_nonempty we purposly move through to the
			// intersection union case to see whether the read is a 'union' type of match
		}
	case intersect_union:
		{
			size_t Ngenes = genes.size();

			if (Ngenes == 1)
			{
				//	If we only match to one gene then the answer is simple
				result = &genes.begin()->first;
				type = &genes.begin()->second.begin()->first;
			}
			else if (Ngenes > 1)
			{
				if (countMode == intersect_union)
				{
					//	For union mode, matching to multiple genes is classed as amiguous
					result = &ambiguousString;
				}
				else 
				{
					//  but for intersection nonempty we can ignore the segments that match nothing and look at all the genes
					//	which match all of the rest of the remaining segments.  We select the gene where the 
					//	length of the match is longest
					size_t bestLength = 0;

					result = &noFeatureString;

					for (auto & gene : genes)
					{
						for (auto & regionType : gene.second)
						{
							if (regionType.second.partial == (genes.nSegments - genes.noMatch))
							{
								if (regionType.second.length > bestLength)
								{
									result = &gene.first;
									type = &regionType.first;

									bestLength = regionType.second.length;
								}
								else if (regionType.second.length == bestLength)
								{
									//	Two genes with the same match length
									result = &ambiguousString;
									type = &blankString;
								}
							}
						}
					}
					//	Spcial case where one of the segments overlap two genes, but the overlaps do not 
					//	fully overlap each other.  This is treated as no feature, and not ambiguous
					if (nonOverlappingGenes && (result != &ambiguousString))
					{
						result = &noFeatureString;
						type = &blankString;
					}

					mode = &nonemptyString;
					break;
				}
			}
			else
			{
				result = &noFeatureString;
			}

			mode = &unionString;
			break;
		}
	}

	if (outputFile.is_open())
			outputFile.printEnd(*mode,*result,*type,segments.name);
	geneCounts[*result][*type]++;

}


void LiBiCount::incBamCounter(const BamAlignment * ba,size_t size)
{
	if ((++bamCounter % 100000) == 0)
		if (verbose)
		{
			cerr << bamCounter << " BAM alignment record pairs processed.";
			if (ba)
				cerr << " cache size = " << size << "  " << references[ba->RefID].RefName << ":" << ba ->Position << endl;
			else
				cerr << " cache reads processed = " << size << endl;
		}
}

bool LiBiCount::processOrderedBamData()
{
	BamAlignment ba1;
	BamAlignment ba2;

	size_t misPairs(0);
	bamCounter = 0;

	bool OK = reader.GetNextAlignment(ba1,false);

	while (OK)
	{
		_DBG(string name = ba1.Name;
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1102:9098:8124");)

			
		regionLists regions;
		regions.name = ba1.Name;

		bool readAlreadyRead = false;
		
		regions.GetRegions(ba1);

		if (ba1.IsFirstMate())
		{
			reader.GetNextAlignment(ba2,false);

			if (ba2.Name == ba1.Name)
			{
				regions.GetRegions(ba2);
			}
			else
			{
				if (bamCounter < 100)
				{
					if (ba1.IsMateMapped())
					{
						if (misPairs++ > 20)
						{
							if (verbose)
								cerr << misPairs << " missing pairs found in the first " << bamCounter << " reads, so assuming data is not name ordered" << endl;
							return false;
						}
					}
				}
				readAlreadyRead = true;
			}
		}

		incBamCounter();

		addRead(regions,genomeDef);

		if (readAlreadyRead)
			swap(ba1,ba2);
		else
			OK = reader.GetNextAlignment(ba1);
	}
	return true;
}


bool LiBiCount::processUnorderedBamData()
{
	BamAlignment ba;

	bamCounter = 0;
	int cacheCounter = 0;
	class readCache : public map<string,regionLists>
	{
	public:
		void save(const string & filename)
		{
			TsvFile outFile;
			outFile.open(filename);
			for (auto & i : This)
				outFile.print(i.first,i.second);
			clear();
		}
	} readCache;

	bool OK = reader.GetNextAlignment(ba,false);


	while (OK)
	{
		_DBG(string name = ba.Name;
		bool found = (name == "HWI-D00133:18:DTWTJACXX:4:1103:6352:25943");)

		bool readAlreadyRead = false;

		if (ba.IsMapped())
		{
			if (ba.IsMateMapped())
			{
				map<string,regionLists>::iterator i = readCache.find(ba.Name);
				if (i == readCache.end())
				{
					regionLists regions;
					regions.GetRegions(ba);
					readCache.emplace(ba.Name,move(regions));
				}
				else
				{
					regionLists & regions = i->second;
					regions.name = ba.Name;

					regions.GetRegions(ba);

					incBamCounter(&ba,readCache.size());

					addRead(regions,genomeDef);

					readCache.erase(i);
				}
			}
			else
			{
				regionLists regions;
				regions.name = ba.Name;
				regions.GetRegions(ba);

				incBamCounter(&ba,readCache.size());

				addRead(regions,genomeDef);

			}
		}
		else if (!ba.IsPaired() || (!ba.IsMateMapped() && ba.IsFirstMate()))
		{
			incBamCounter(&ba,readCache.size());
		}
		if (readCache.size() > maxCacheSize)
		{
			if (verbose)
				cerr << "Outputting cache data " << cacheCounter+1 << endl;
			readCache.save(resultsFilename.replaceSuffix(".temp.",cacheCounter++));
		}

		OK = reader.GetNextAlignment(ba,false);
	}

	if (cacheCounter == 0)
	{
		size_t cacheReadCounts = 0;
		for (auto & i : readCache)
		{
			i.second.name = i.first;
			addRead(i.second,genomeDef);
			incBamCounter(0,cacheReadCounts++);
		}
	}
	else
	{
		readCache.save(resultsFilename.replaceSuffix(".temp.",cacheCounter++));
		processCachedReads(cacheCounter);
	}

	return true;
}

void LiBiCount::processCachedReads(size_t cacheFileCount)
{
/*	Go through the reads in the file caches.  The read pairs associated with the same fragment will be in different files
	and the files are ordered by read name so we only need to look at the next reads in each file to spot pairs that can be processed
	as a pair
	*/

	vector<cacheRead> cacheReads(cacheFileCount);

	//	readIndex has a lits of the current reads, ordered by name.
	multimap<stringEx,int> readIndex;

	for (int i = 0;i < cacheFileCount;i++)
	{
		cacheReads[i].open(resultsFilename.replaceSuffix(".temp.",i));
		readIndex.emplace(cacheReads[i].name,i);
	}

	int cacheReadCounter = 0;
	while (readIndex.size())
	{
		_DBG( bool found = (readIndex.begin()->first == "HWI-D00133:18:DTWTJACXX:4:1103:6352:25943");)

		auto i1 = readIndex.begin();
		auto i2 = next(i1,1);

		int index1 = i1->second;

		if ((i2 != readIndex.end()) && (i1->first == i2->first))
		{
			//	We have the two ends of a paired end read.  Combine them and calculate counts
			int index2 = i2->second;
			cacheReads[index1].combine(cacheReads[index2]);
			if (i1->first)								//Avoid adding null entries with no read name
				addRead(cacheReads[index1],genomeDef);
			//	Erase the read and get the next one from the associated cache files
			readIndex.erase(i2);
			if(cacheReads[index2].readNext())
				readIndex.emplace(cacheReads[index2].name,index2);
		}
		else
		{
			if (i1->first)
				addRead(cacheReads[index1],genomeDef);
		}

		readIndex.erase(i1);
		if(cacheReads[index1].readNext())
			readIndex.emplace(cacheReads[index1].name,index1);
		incBamCounter(0,++cacheReadCounter);

	}
	cacheReads.clear();

	for (int i = 0;i < cacheFileCount;i++)
	{
		remove(resultsFilename.replaceSuffix(".temp.",i).c_str());
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



cacheRead::~cacheRead() 
{
	close();
};


bool cacheRead::open(const std::string filename)
{
	file = new std::ifstream();
	file ->open(filename);
	if (!file ->is_open()) return false;
	readNext();
	return true;
}
bool cacheRead::readNext()
{
	if (file->eof())
		return false;
	data.clear();
	std::string line;
	getline(*file,line);
	parseTsv(line,name,data);
	return true;
}
void cacheRead::close()
{
	if (file)
	{
		file->close();
		delete (file);
		file = 0;
	}
}
