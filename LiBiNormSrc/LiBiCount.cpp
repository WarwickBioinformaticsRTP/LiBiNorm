#include <stdlib.h>

#ifdef _WIN32
#include <crtdbg.h>
//	for rmdir and mkdir
#include <direct.h>
#else
//	For rmdir
#include <unistd.h>	
//	for mkdir
#include <sys/stat.h>
#define mkdir(A) mkdir(A,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

#include "libCommon.h"
#include "LiBiCount.h"
#include "parser.h"

#define MATCH_USING_POSITION

#ifdef _DEBUG
#define READ_CACHE_SIZE 50
#define REP_LEN 100
#else
#define READ_CACHE_SIZE 2000000
//#define READ_CACHE_SIZE 10000
#define REP_LEN 100000
#endif

using namespace std;

#define BAMNAME "HISEQ2000-05:531:C7LLMACXX:2:1101:16672:4313"

int LiBiCount::main(int argc, char **argv)
{

	stringEx bamFileName,gtfFileName,outputFilename,
		id_attribute = "gene_id";

	setEx<string> feature_type;
	bool featureAdded = false;

	reverseStrand = false;
	useStrand = true;
	verbose = true;
	htSeqCompatible = true;
	minqual = 10;
	nameOrder = true;
	countMode = intersect_union;
	maxCacheSize = READ_CACHE_SIZE;

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1],"-h") ==0) || (strcmp(argv[1],"--help")==0))))
	{
printf("Usage: LiBiNorm [options] alignment_file gff_file\n");
printf(" based on htseq-count\n");
printf("This program takes an alignment file in SAM/BAM format and a feature file in\n");
printf("GFF format and calculates for each feature the number of reads mapping to it.\n");
printf("See http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details.\n");
printf("\n");
printf("Options:\n");
printf("  -h, --help            show this help message and exit\n");
printf("  -r ORDER, --order=ORDER\n");
printf("                        'pos' or 'name'. Sorting order of <alignment_file>\n");
printf("                        (default: name). Paired-end sequencing data must be\n");
printf("                        sorted either by position or by read name, and the\n");
printf("                        sorting order must be specified. Ignored for single-\n");
printf("                        end data.\n");
printf("  -s STRANDED, --stranded=STRANDED\n");
printf("                        whether the data is from a strand-specific assay.\n");
printf("                        Specify 'yes', 'no', or 'reverse' (default: yes).\n");
printf("                        'reverse' means 'yes' with reversed strand\n");
printf("                        interpretation\n");
printf("  -a MINAQUAL, --minaqual=MINAQUAL\n");
printf("                        skip all reads with alignment quality lower than the\n");
printf("                        given minimum value (default: 10)\n");
printf("  -t FEATURETYPE, --type=FEATURETYPE\n");
printf("                        feature type (3rd column in GFF file) to be used, all\n");
printf("                        features of other type are ignored (default, suitable\n");
printf("                        for Ensembl GTF files: exon)\n");
printf("  -i IDATTR, --idattr=IDATTR\n");
printf("                        GFF attribute to be used as feature ID (default,\n");
printf("                        suitable for Ensembl GTF files: gene_id)\n");
printf("  -m MODE, --mode=MODE  mode to handle reads overlapping more than one feature\n");
printf("                        (choices: union, intersection-strict, intersection-\n");
printf("                        nonempty; default: union)\n");
printf("  -c COUNTS, --counts=COUNT\n");
printf("                        Name of output file. default: writes to stdout)\n");
//printf("  -o SAMOUT, --samout=SAMOUT\n");
//printf("                        write out all SAM alignment records into an output SAM\n");
//printf("                        file called SAMOUT, annotating each line with its\n");
//printf("                        feature assignment (as an optional field with tag\n");
//printf("                        'XF')\n");
printf("  -q, --quiet           suppress progress report\n");
printf("\n");
printf("Written by Nigel Dyer (nigel.dyer@warwick.ac.uk)\n");
		return EXIT_SUCCESS;
	}

	int ni = 1;

	if (argc < 3)
		exitFail("Insufficient arguments");

	while(ni < argc-2)
	{
		bool opt2 = false;
/*		if(strcmp(argv[ni], "-t") == 0)
		{
			fileCompare(argc - 1,&argv[1]);
			exitSuccess();
		}
		else*/ if((strcmp(argv[ni], "-h") == 0) || (opt2 = (strncmp(argv[ni], "--cache=",8) == 0)))
		{
			maxCacheSize = atoi(opt2?argv[ni]+8:argv[++ni]);
		}
		else if ((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--stranded=",11) == 0)))
		{
			string strand(opt2?argv[ni]+11:argv[++ni]);
			if (strand == "yes")
			{}
			else if (strand == "reverse")
				reverseStrand = true;
			else if (strand == "no")
				useStrand = false;
			else
				exitFail("Invalid strand option: ",strand);
		}
		else if ((strcmp(argv[ni], "-a") == 0) || (opt2=(strncmp(argv[ni], "--minaqual=",11) == 0)))
		{
			minqual = atoi(opt2?argv[ni]+11:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-f") == 0) || (opt2=(strncmp(argv[ni], "--type=",7) == 0)))
		{
			feature_type.emplace(opt2?argv[ni]+7:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-i") == 0) || (opt2=(strncmp(argv[ni], "--idattr=",9) == 0))) 
		{
			id_attribute = opt2?argv[ni]+9:argv[++ni];
		}
		else if((strcmp(argv[ni], "-r") == 0)|| (opt2=(strncmp(argv[ni], "--order=",8) == 0)))
		{
			string mode(opt2?argv[ni]+8:argv[++ni]);
			if (mode == "pos")
				nameOrder = false;
			else if (mode == "name")
				nameOrder = true;
			else
				exitFail("Unknown 'order':",mode,".  Should be pos or name");
		}
		else if ((strcmp(argv[ni], "-m") == 0) || (opt2 = (strncmp(argv[ni], "--mode=",7) == 0)))
		{
			string mode = opt2?argv[ni]+7:argv[++ni];
			if (mode == "union")
				countMode = intersect_union;
			else if (mode == "intersection-strict")
				countMode = intersect_strict;
			else if (mode == "intersection-nonempty")
				countMode = intersect_nonempty;
			else if (mode == "intersection-all")
				countMode = intersect_all;
			else
				exitFail("Invalid mode: ",mode);
		}
		else if((strcmp(argv[ni], "-q") == 0) || (strcmp(argv[ni], "--quiet") == 0))
		{
			verbose = false;
		}
/*		else if((strcmp(argv[ni], "-o") == 0) || (strcmp(argv[ni], "--samout") == 0))
		{
			outputFilename = argv[++ni];
		}*/
		else if((strcmp(argv[ni], "-c") == 0) || (opt2 = (strncmp(argv[ni], "--counts=",9) == 0)))
		{
			countsFilename = opt2?argv[++ni]+9:argv[++ni];
			tempDirectory = countsFilename.replaceSuffix("_tempFiles");
		}
		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	bamFileName = argv[argc-2];
	gtfFileName = argv[argc-1];


	if (feature_type.size() == 0)
		feature_type.emplace("exon");

	if (outputFilename && !outputFile.open(outputFilename))
		exitFail("Unable to open output file: ",outputFilename);

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	if (!tempDirectory)
	{
		const char * p = getenv("TEMP");
		if (p == 0) 
		   p = getenv("TMPDIR");

		if (p == 0)
			exitFail("Unable to identify temporary directory from environment variable, Try using the -c option for the count files instead");  	

		tempDirectory = p;
		tempDirectory += stringEx("/LiBiNorm_temp_",rand(),rand());

		if (verbose)
			cerr << "temp Directory = " << tempDirectory << endl;
	}

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();

	initClock();

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

	if (nameOrder)
		processNameOrderedBamData();
	else
	{
		//Only need temporary directory for caching files
		// if the data are position ordered
		mkdir(tempDirectory.c_str());
		tempDirectory += "/";
		processPositionOrderedBamData();
	}

	if(!outputGeneCounts(countsFilename))
		exitFail("Unable to output counts to :",countsFilename);

	if (verbose)
		elapsedTime();

	string test;
	_DBG(cin >> test);
//	cin >> test;

	if (!nameOrder)
		rmdir(tempDirectory.c_str());


	return EXIT_SUCCESS;
}



bool LiBiCount::outputGeneCounts(const string & filename)
{
	TsvFile output;
	
	if (!output.open(filename))
		return false;

	for(auto i : geneCounts)
	{
		if (i.first.substr(0,2) != "__")
			geneCounts.print(i.first,output);
/*		for(auto j : i.second)
		{
			if (i.first.substr(0,2) != "__")
				output.printEnd(i.first,j.first,j.second);
		}*/
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
	if (segments.NH > 1)
	{
		geneCounts[notUnique]++;
		return;
	}
	else if (segments.qual < minqual)
	{
		geneCounts[lowQualString]++;
		return;
	}


	//	The paired end read consists of a number of segments.   If the two ends were aligned to different chromosomes
	//	then the segments will be on different chromosomes
	chromosomeGeneInfo genes;

	const string * chromosome = &blankString;
	string location;

	size_t pairNo = 0;
	for (auto & chromSegments : segments.data)
	{
		//	For teh segments on each of the chromosomes (normally only one chromosome) get the gtfRegions for the chromosome
		genomeGtfRegions::const_iterator thisChromGtfRegions = gtfData.genomeGtfData.find(references[chromSegments.first].RefName);

		if (thisChromGtfRegions != gtfData.genomeGtfData.end())
		{
			chromosome = &references[chromSegments.first].RefName;

			if (outputFile.is_open() && location.empty())
			{
				location = stringEx(*chromosome,":",chromSegments.second.data.begin()->first);
			}

			//	Get the map of neds of gtfRegions associated with the chromosome
			const chromosomeEndIndexMap & thisChromEndMap = gtfData.genomeEndIndex.at(references[chromSegments.first].RefName);

			//	And find the first one that finishes at or beyond the start of the first segment
			chromosomeEndIndexMap::const_iterator indirectIteratorStart = thisChromEndMap.lower_bound(chromSegments.second.data.begin()->second.start);

			//	And move back one to ensure we have the region that covers segment
//			if (indirectIteratorStart != thisChromEndMap.begin())
			if (indirectIteratorStart == thisChromEndMap.end())
				indirectIteratorStart--;

			chromosomeGtfData::iterator gtfRegion = indirectIteratorStart->second;

			//	If this region overlaps any other regions then go to the one that starts the earliest.
			//	If there were no overlaps then default is the overlaps points to self
			gtfRegion = *gtfRegion->second.overlaps;


			//And now go through each of the segments
			for (auto & segment : chromSegments.second.data)
			{

				genes.nSegments++;
				vector<gtfOverlap> overlaps;
				//	Trying out each of the gtfRegions in turn to see if there is an overlap. 
				//	If there is then it gets added to the list of ovberlaps
				for (auto j = gtfRegion; (j != thisChromGtfRegions->second.end()) && (j->first <= segment.second.end); j++)
				{
					//	Check for strand match
					if (!useStrand || ((j->second.strand == segment.second.strand) != reverseStrand))
					{
						j->second.checkOverlap(segment.second,overlaps);
					}

					gtfRegion = *j->second.overlaps;
				}


				//	Now go through all the overlaps between the read and the regions identified in the gtf file
				if (overlaps.size() == 0)
				{
					//	There were none
					genes.noMatch++;
				}
				else
				{
					struct gtfId
					{
						string name,type;
						gtfId(){};
						gtfId(const gtfOverlap & region):name(region.geneName),type(region.featType){};

					};

					class segRegion
					{
					public:
						size_t min,max;
						vector<gtfId> geneSet;
						segRegion():min(INT_MAX),max(0){};
						segRegion(const gtfOverlap & region):min(region.start),max(region.finish) 
						{
							geneSet.emplace_back(region);
						};

					};

					//	
					gtfOverlap & overlap = overlaps[0];

					//	This will create an entry for the combination if it does not exist before, which is needed later on
					overlapCounts & GeneAttributeCombo1 = genes[overlap.geneName][overlap.featType];

					//	Need to register all strict overlaps, because overlaps-strict require them to be unique
					if (overlap.strict)
						GeneAttributeCombo1.strict++;

					if (overlaps.size() == 1)
					{
						GeneAttributeCombo1.length += (overlap.finish - overlap.start);
						GeneAttributeCombo1.partial++;
					}
					else
					{
						vector<segRegion> segRegions;
						segRegions.emplace_back(overlap);

						// More than one region.  If the gtf regions identify separate sections of the read then keep both, 
						//	If they overlap then only keep the largest if it fully overlaps the other.  If they are identical then keep both
						//	

						for (int i = 1; i < overlaps.size(); i++)
						{
							gtfOverlap & overlap = overlaps[i];

							//	Create dummy entry for every region type that the read overlapped;
		
							overlapCounts & GeneAttributeCombo2 = genes[overlap.geneName][overlap.featType];

							if (overlap.strict)
								GeneAttributeCombo2.strict++;
							bool newRegion = false;

							for (auto & region : segRegions)
							{
								if ((overlap.start == region.min) && (overlap.finish == region.max))
								{
									//Identical
									region.geneSet.emplace_back(overlap);
									newRegion = false;
								}
								else if ((overlap.start >= region.min) && (overlap.finish <= region.max))
								{
									//	smaller ignore.  It must be smaller in that we have already excluded the case
									//	where it is identical
									newRegion = false;
								}
								else if ((overlap.start <= region.min) && (overlap.finish >= region.max))
								{
									//	It is bigger so replace.  Again we have excluded the option that it is identical
									//	which would have been otherwise included in the case
									region.min = overlap.start;
									region.max = overlap.finish;
									region.geneSet.resize(1);
									region.geneSet.at(0).name = overlap.geneName;
									region.geneSet.at(0).type = overlap.featType;
									newRegion = false;
								}
								else
								{
									newRegion = true;
								}
							}
							if (newRegion)
							{
								segRegions.emplace_back(overlap);
							}
						}

						//	We now have one or more regions within the read, each one of which matches regions in the gtf file
						//	
						for (gtfId & g : segRegions.at(0).geneSet)
						{
							bool matchesInAllRegions = true;
							int size = segRegions.at(0).max - segRegions.at(0).min;
							for (size_t i = 1; i < segRegions.size();i++)
							{
								bool found = false;
								for (auto & j: segRegions.at(i).geneSet)
								{
									if ((g.name == j.name) && (g.type == j.type))
									{
										size += (segRegions.at(i).max - segRegions.at(i).min);
										found = true;
										break;
									}
								}
								if (!found)
									matchesInAllRegions = false;
							}
							if (matchesInAllRegions)
							{
								genes.at(g.name).at(g.type).partial++;
								genes.at(g.name).at(g.type).length += size;
							}
						}
					}
				}
			}
		}
		else
		{
			//	There were no genes identified on this 'chromosome' of the reference sequence
			if (htSeqCompatible)
			{
				genes.clear();
			}
			else
			{
				genes.nSegments++;
				genes.noMatch++;
			}
		}
		pairNo++;
	}


	//	Result mode options
	static string strictString = "strict";
	static string nonemptyString = "nonempty";
	static string unionString = "union";

	//	Use pointers to strings rather than the strings themselves for efficiency as the it avoids
	// creating and deleting copies of strings
	const string * result = nullptr;
	const string * type = &blankString;		//
	const string * mode = &blankString;		//The mode that selected the region

	switch(countMode)
	{
	case intersect_strict:
	case intersect_nonempty:		
	case intersect_all:		
		{
			mode = &strictString;

			//	For a strict match, all of the read segments must lie inside an annotated region of the same gene
			for (auto & gene : genes)
			{
				for (auto & regionType : gene.second)
				{
					if (regionType.second.strict == genes.nSegments)
					{
						if (result)
						{
							if (countMode == intersect_all)
							{
								//	In intersect all we include all of the options, ie there will be multiple counts associated with
								//	a single fragment
								if (outputFile.is_open())
									outputFile.printEnd(*mode,*result,*type,location,segments.name);
								geneCounts.at(*result).at(*type)++;

								result = &gene.first;
								type = &regionType.first;
							}
							else
							{
								//	If we have two strict matches then the result is ambigous, no need to look any further
								result = &ambiguousString;
								type = &blankString;
								break;
							}
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

			if (result == nullptr)
			{
				if ((countMode == intersect_strict) || (countMode == intersect_all))
				{
					result = &noFeatureString;
					break;
				}
				//	otherwise it is nonempty mode, go on to see if there is a union style match
			}
			else
			{
				break;
			}

			//   ***************   BEWARE******************
			//	If the mode is intersection_nonempty we purposly move through to the
			// intersection union case to see whether the read is a 'union' type of match
			//	so, no break between case statements
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
		outputFile.printEnd(*mode,*result,*type,location,segments.name);

	geneCounts.at(*result).at(*type)++;

}


void LiBiCount::incBamCounter(const BamAlignment * ba,size_t size)
{
	if ((++bamCounter % REP_LEN) == 0)
		if (verbose)
		{
			cerr << bamCounter << " BAM alignment record pairs processed.";
			if (ba)
			{
				if(ba->RefID >= 0)
					cerr << " cache size = " << size << "  " << references[ba->RefID].RefName << ":" << ba ->Position << endl;
				else
					cerr << " cache size = " << size << "  unmapped read" << endl;
			}
			else if (size != -1)
				cerr << " cache reads processed = " << size << endl;
			else
				cerr << endl;
		}
}

bool LiBiCount::AReadIsMapped(const BamAlignment & ba)
{

	if (!ba.IsMapped() && !ba.IsMateMapped())
	{
		if ((ba.IsPaired() && ba.IsFirstMate()) || !ba.IsPaired())
			geneCounts[notAlignedString]++;
		return false;
	}

	return true;
}

#define READ_BUFFER_SIZE 100

bool LiBiCount::processNameOrderedBamData()
{
	BamAlignment ba[READ_BUFFER_SIZE];
	bool used[READ_BUFFER_SIZE];

	size_t misPairs(0);
	bamCounter = 0;
	set<string> previousNames;

	bool OK = reader.GetNextAlignment(ba[0],false);

	int N = 0;

	while (OK)
	{
		string & name = ba[0].Name;
		_DBG(bool found = (name == BAMNAME);)

		int Nreads = 0;

		do {
			used[Nreads] = false;
			OK = reader.GetNextAlignment(ba[++Nreads],false);
		}
		while ((ba[Nreads].Name == name) && (OK) && (Nreads < READ_BUFFER_SIZE));

		//	This is a fairly direct implementation of the htseq logic in __init_.py line 570 onwards
		//	There are potential issues with the ability of teh code to cope with datasets where there are multiple
		//	alignments

		if (AReadIsMapped(ba[0]))
		{
			for (size_t i = 0;i < Nreads;i++)
			{
				if (!used[i])
				{
					regionLists regions(readData(move(ba[i])),name);
					for (size_t j = i+1;j < Nreads;j++)
					{
						if (!used[j])
						{
							if ((ba[i].IsFirstMate() == ba[j].IsFirstMate())
								|| (ba[i].IsMapped() != ba[j].IsMateMapped())
								|| (ba[i].IsMateMapped() != ba[j].IsMapped()))
								continue;

							if (!(ba[i].IsMapped() && ba[j].IsMapped())  ||
									
								((ba[i].RefID == ba[j].MateRefID)
								&& (ba[i].MateRefID == ba[j].RefID)
								&& (ba[i].Position == ba[j].MatePosition)
								&& (ba[i].MatePosition == ba[j].Position)
								))
							{
								used[i] = true;
								used[j] = true;		
								regions.combine(move(ba[j]));
								addRead(regions,genomeDef);
								incBamCounter();
								break;
							}
						}
					}
					if (!used[i])
					{
						addRead(regions,genomeDef);
						incBamCounter();
					}
				}
			}
		}

		swap(ba[0],ba[Nreads]);

	}
	return true;
}


bool LiBiCount::processPositionOrderedBamData()
{
	BamAlignment ba;

	bamCounter = 0;
	int cacheCounter = 0;


	class readCacheClass : public map<string,vector<readData> >
	{
	public:
		void save(const string & filename)
		{
			TsvFile outFile;
			outFile.open(filename);
			for (auto & i : This)
				for (auto & j : i.second)
					outFile.print(i.first,j);
			clear();
		}
	} readCache;

	bool OK = reader.GetNextAlignment(ba,false);

	while (OK)
	{
		_DBG(string name = ba.Name;
		bool found = (name == BAMNAME);)

		//	The NH handling is complex is that there may be one NH (with NH = 1) at one end, and multiple NHs (with NH > 1) at the other
		//	The NH > 1 samples have to be used to either pair with the other, or to remove the NH = 1 sample 

		if (AReadIsMapped(ba))
		{
			//	Store reads in a cache so they can be paired up.
			stringEx mateIndex(ba.Name,"_",
#ifdef MATCH_USING_POSITION
				ba.IsMateMapped()?stringEx(ba.MateRefID,ba.MatePosition):"0","_",
				ba.IsMapped()?stringEx(ba.RefID,ba.Position):"0","_",
#endif
				(ba.IsMapped() && ba.IsMateMapped())?-ba.InsertSize:0,"_",
				ba.IsFirstMate()?"S":"F");
			readCacheClass::iterator i = readCache.find(mateIndex);
			if (i == readCache.end())
			{
				stringEx thisIndex(ba.Name,"_",
#ifdef MATCH_USING_POSITION
					ba.IsMapped()?stringEx(ba.RefID,ba.Position):"0","_",
					ba.IsMateMapped()?stringEx(ba.MateRefID,ba.MatePosition):"0","_",
#endif
					(ba.IsMapped() && ba.IsMateMapped())?ba.InsertSize:0,"_",
					ba.IsFirstMate()?"F":"S");

				i = readCache.find(thisIndex);
				if (i == readCache.end())
				{
					auto i = readCache.emplace(thisIndex,vector<readData>());
					i.first->second.emplace_back(move(ba));
				}
				else
					i->second.emplace_back(move(ba));
			}
			else
			{
				regionLists regions(i->second.front(),ba.Name);

				regions.combine(move(ba));

				addRead(regions,genomeDef);

				if (i->second.size() > 1)
					i->second.erase(i->second.begin());
				else
					readCache.erase(i);

				incBamCounter(&ba,readCache.size());

			}
		}
		incBamCounter(&ba,readCache.size());

		if (readCache.size() > maxCacheSize)
		{
			if (verbose)
				cerr << "Outputting cache data " << cacheCounter+1 << endl;
			readCache.save(stringEx(tempDirectory,"file",cacheCounter++));
		}

		OK = reader.GetNextAlignment(ba,false);
	}

	if (cacheCounter == 0)
	{
		//	There are no records cached to disk, so the internal cache simply contains unpaired reads that can be processed
		//	immediately
		size_t cacheReadCounts = 0;
		for (auto & i : readCache)
		{
			_DBG(
				vector<string> tags;
				parser(i.first,"_",tags);
				bool found = (tags[0] == BAMNAME);
				)

			for (auto & j: i.second)
			{

				regionLists rl(j,i.first);
				addRead(rl,genomeDef);

				incBamCounter(0,cacheReadCounts++);
			}
		}
	}
	else
	{
		//	We have cached some records to disk, so the internal cache must be flushed as well so that all of the cached records can be
		//	processed as an ensemble
		readCache.save(stringEx(tempDirectory,"file",cacheCounter++));
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

	vector<cacheEntry> cacheReads(cacheFileCount);

	//	readIndex has a lits of the current reads, ordered by name.

	class readCache : public map<string,map<int,vector<readData> > > 
	{
	public:
		void replace(iterator & i,vector<cacheEntry> & cacheReads)
		{
			int index = i->second.begin()->first;
			vector<readData> & r = i->second.begin()->second;
			if (r.size() > 1)
				r.erase(r.begin());
			else
			{
				map<int,vector<readData> > & s = i->second;
				if (s.size() > 1)
					s.erase(s.begin());
				else
					erase(i);
			}

			if(cacheReads[index].readNext())
				This[cacheReads[index].name][index].emplace_back(cacheReads[index]);

		}
	} reads;

	for (int i = 0;i < cacheFileCount;i++)
	{
		cacheReads[i].open(stringEx(tempDirectory,"file",i));
		reads[cacheReads[i].name][i].emplace_back(cacheReads[i]);
		for (int j = 0;(j < 10) && (cacheReads[i].readNext());j++)
		{
			reads[cacheReads[i].name][i].emplace_back(cacheReads[i]);
		}
	}

	int cacheReadCounter = 0;
	while (reads.size())
	{
		_DBG( bool found = (reads.begin()->first == BAMNAME);)

		readCache::iterator i1 = reads.begin();

		size_t nameLen = i1->first.length();
		if (nameLen == 0)
		{
			reads.erase(i1);
		}
		else
		{
			string name = i1->first;
			vector<string> nameParts;

			parser(name,"_",nameParts);

			name = stringEx(nameParts[0],
#ifdef MATCH_USING_POSITION
				"_",nameParts[2],"_",nameParts[1],"_",-atoi(nameParts[3].c_str()),"_",(nameParts[4] == "F")?"S":"F");
#else
				"_",-atoi(nameParts[1].c_str()),"_",(nameParts[2] == "F")?"S":"F");
#endif

			regionLists rl(i1->second.begin()->second[0],nameParts[0]);

			readCache::iterator i2 = reads.find(name);
			if (i2 != reads.end())
			{
				//	We have the two ends of a paired end read.  Combine them and calculate counts
				rl.combine(i2->second.begin()->second[0]);
				reads.replace(i2,cacheReads);
			}

			addRead(rl,genomeDef);
			reads.replace(i1,cacheReads);

			incBamCounter(0,++cacheReadCounter);
		}
	}
	
	//	Closes all of the files and then deletes them
	cacheReads.clear();
}


void LiBiCount::fileCompare(int argc, char **argv)
{
	string matches("Y:\\SysmedIBD\\CD\\test\\matches.sam"),
		myMatches("Y:\\SysmedIBD\\CD\\test\\mymatches.txt"),
		comparison("Y:\\SysmedIBD\\CD\\test\\comparison.txt");

	int ni = 1;
	while(ni < argc)
	{
		switch (ni)
		{
		case 1: matches = argv[ni++]; break;
		case 2: myMatches = argv[ni++]; break;
		case 3: comparison = argv[ni++]; break;
		}
	}



	ifstream samFile;
	samFile.open(matches);
	if (!samFile.is_open())
		exitFail("Failed to open samfile: ",matches);

	ifstream myFile;
	myFile.open(myMatches);
	if (!myFile.is_open())
		exitFail("Failed to open LiBiNorm output file: ",myMatches);



	TsvFile outFile;
	if (!outFile.open(comparison))
		exitFail("Failed to open output file: ",comparison);

	string line;

	vectorEx<size_t> cols(1,2,3,4,5);


	vector<string> myParams;
	vector<stringEx> samParams[2];

	bool OK = getline(samFile,line).good();

	parseTsv(line,samParams[0]);
	int Entry = 0;

	while (OK) {

		Entry++;
		getline(myFile,line);
		myParams.clear();
		parseTsv(line,myParams);


		getline(samFile,line);
		samParams[1].clear();
		parseTsv(line,samParams[1]);

		if (samParams[0][0] != myParams[4])
			exitFail("Mismatch entry ",Entry,"  ",samParams[0][0],"\n",myParams[4]);

		string match;
		for (size_t i = 8;i < samParams[0].size();i++)
		{
			if (samParams[0][i].startsWith("XF:Z:"))
			{
				match = samParams[0][i].substr(5);
			}
		}


		if (!(((match.substr(0,2) == "__") && (myParams[1].substr(0,2) == "__")) || (match == myParams[1])))
		{
			outFile.printStart(myParams[4],myParams[0],myParams[1],myParams[2]);
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
		}

		if (samParams[1][0] != myParams[4])
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



cacheEntry::~cacheEntry() 
{
	close();
};


bool cacheEntry::open(const std::string filename)
{
	file = new std::ifstream();
	file ->open(filename);
	if (!file ->is_open()) return false;
	fname = filename;
	readNext();
	return true;
}
bool cacheEntry::readNext()
{
	if (file->eof())
		return false;
	cigar.clear();
	std::string line;
	getline(*file,line);
	parseTsv(line,name,refId,position,strand,cigar,NH,qual);
	return true;
}
void cacheEntry::close()
{
	if (file)
	{
		file->close();
		delete (file);
		file = 0;
		remove(fname.c_str());
	}
}
