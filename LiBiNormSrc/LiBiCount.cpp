

#ifdef _WIN32
#include <direct.h>
#else
//	For rmdir
#include <unistd.h>	
#endif

#include "libCommon.h"
#include "LiBiCount.h"
#include "LiBiNorm.h"
#include "parser.h"

//#define MATCH_USING_POSITION

//	The insert size for matching pairs should be A and -A.  In some datasets they are A and A.  By using the absolute value of the
//	insert size we can ensure that the pairs are still matched up
#define USE_ABS_INSERT_TO_MATCH_READS


#ifdef USE_ABS_INSERT_TO_MATCH_READS
#define insertConv(A) abs(A)
#else
#define insertConv(A) A
#endif


#ifdef _DEBUG
//	Put reads into cache file when number of reads exceed READ_CACHE_SIZE
#define READ_CACHE_SIZE 50000
//	Report progress every REP_LEN entries
#define REP_LEN 100000
#define BAMNAME "SRR557798.3975438"
bool dbgFound = false;
#else
#define READ_CACHE_SIZE 2000000
//#define READ_CACHE_SIZE 10000
#define REP_LEN 100000
#endif


#ifndef _DEBUG

//	Test code: Compare lengths calculated from the gff file with the lengths in the original landscape file
//#define COMPARE_RESULTS "Y:\\LiBiNorm\\SRR557798\\SRR557798.NoA.plus.minus"

//  Output detailed results of interpreting the Feature file
//#define OUTPUT_FEATURE_DATA
#endif

using namespace std;


int LiBiCount::main(int argc, char **argv)
{

	stringEx bamFileName,featureFileName,outputFilename,geneListFilename,
		id_attribute = "gene_id";

	stringEx feature_type;

	reverseStrand = false;
	useStrand = true;
	verbose = true;
	minqual = 10;
	nameOrder = true;
	countMode = DEFAULT_COUNT_MODE;
	maxCacheSize = READ_CACHE_SIZE;

	if(argc < 1)
	{
		exitFail("Error: parameter wrong!");
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
printf("                        features of other type are ignored (default for \n");
printf(_s("                        Ensemble GTF and GFF files: ", DEFAULT_FEATURE_TYPE_EXON,")\n"));
printf("  -i IDATTR, --idattr=IDATTR\n");
printf("                        GFF attribute to be used as feature ID (default for \n");
printf(_s("                        Ensembl GTF files: ", DEFAULT_GTF_ID_ATTRIBUTE, "\n"));
printf(_s("                        default for GFF3 files: ", DEFAULT_GFF_ID_ATTRIBUTE, ")\n"));
printf("  -m MODE, --mode=MODE  mode to handle reads overlapping more than one feature\n");
printf("                        (choices: union, intersection-strict, intersection-\n");
printf("                        nonempty; default: union)\n");
printf("  -z, --htseq-compatible\n");
printf("                        Run in htseq-compatible mode\n");

//printf("  -o SAMOUT, --samout=SAMOUT\n");
//printf("                        write out all SAM alignment records into an output SAM\n");
//printf("                        file called SAMOUT, annotating each line with its\n");
//printf("                        feature assignment (as an optional field with tag\n");
//printf("                        'XF')\n");
helpCommon();
printf("  -g F, --genes= F      Only perform the analysis for the genes listed in the file with name F\n");
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
		else if((strcmp(argv[ni], "-k") == 0) || (opt2 = (strncmp(argv[ni], "--cache=",8) == 0)))
		{
			maxCacheSize = atoi(opt2?argv[ni]+8:argv[++ni]);
		}
		else */if ((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--stranded=",11) == 0)))
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
		else if((strcmp(argv[ni], "-t") == 0) || (opt2=(strncmp(argv[ni], "--type=",7) == 0)))
		{
			feature_type = opt2?argv[ni]+7:argv[++ni];
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
		else if ((strcmp(argv[ni], "-z") == 0) || (strcmp(argv[ni], "---htseq-compatible") == 0))
		{
			normalise = false;
			htSeqCompatible = true;
		}
		else if (commandParseCommon(ni, argc-2, argv)) //argc -2 to allow for the two fixed end parameters
		{
		}
		else if ((strcmp(argv[ni], "-g") == 0) || (opt2 = (strncmp(argv[ni], "--genes=", 8) == 0)))
		{
			geneListFilename = opt2 ? argv[++ni] + 8 : argv[++ni];
		}
		else if (strcmp(argv[ni], "-x") == 0)
		{
			pauseAtEnd = true;
		}
		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	bamFileName = argv[argc-2];
	featureFileName = argv[argc-1];

	//	If we have specified -N then we run all of the models for preset number of runs.
	if (normaliseResultsFilename)
	{
		//	Unless specified we find the best model
		NrunsOtherModels = Nruns;
	
	}
	else
	{
		if (theModel == noModel)
			theModel = DEFAULT_MODEL;
	}
	

	if(countsFilename)
		tempDirectory = countsFilename.replaceSuffix("_tempFiles");


	if (!feature_type)
		feature_type = DEFAULT_FEATURE_TYPE_EXON;

	if (!id_attribute)
	{
		if (featureFileName.suffix() == "gtf")
			id_attribute = DEFAULT_GTF_ID_ATTRIBUTE;
		else if (featureFileName.suffix().startsWith("gff"))
			id_attribute = DEFAULT_GFF_ID_ATTRIBUTE;
		else
			exitFail("Unable to identify feature file type in order to specifiy default id attribute");
	}

	if(normalise && (feature_type != "exon"))
		exitFail("Can only normalise data when 'exon' is the only feature specified");

	if (outputFilename && !outputFile.open(outputFilename))
		exitFail("Unable to open output file: ",outputFilename);

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM files: ",bamFileName);

	if (!nameOrder)
		tempDirectory = tempDirectory::get(tempDirectory);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();
	for(auto & i : references)
	{
		if (strncasecmp(i.RefName.c_str(),"chr",3) == 0)
			i.RefName = i.RefName.substr(3);
		genomeDef.addToChromosomeMap(i.RefLength, i.RefName);
	}

	initClock();
	
	if (!genomeDef.open(featureFileName,id_attribute,feature_type))
		exitFail("Could not open feature file: ",featureFileName);

#ifdef	OUTPUT_FEATURE_DATA
		if ((countsFilename) && !genomeDef.printEntries(countsFilename.replaceSuffix("_genome.txt")))
		progMessage("Unable to output genome data to :",countsFilename.replaceSuffix("_genome.txt"));
#endif

	geneCounts.addEntry("reference", DEFAULT_NORMALISATION_GENE_LENGTH);
	if (geneListFilename)
	{
		progMessage("Using genes/transcripts listed in ", geneListFilename);
		geneCounts.useSelectedGenes(geneListFilename);
	}

	genomeDef.index(geneCounts);

#ifdef COMPARE_RESULTS
	transcriptDataMap transData;
	transData.loadData(COMPARE_RESULTS);
	TsvFile testOut;
	testOut.open(countsFilename.replaceSuffix(".compareLengths.txt"));

	for (size_t i = 0; i < transData.size(); i++)
	{
		long len = genomeDef.genes[transData[i].gene].length;
		testOut.print(transData[i].gene, transData[i].length, len);
	}
	testOut.close();

#endif

	elapsedTime("Feature file consolidated.");

	_DBG(genomeDef.outputChromData(featureFileName.replaceSuffix(".txt"));)

	if (nameOrder)
		processNameOrderedBamData();
	else
		processPositionOrderedBamData();

	//	Need to output landscape file now because the data will be modified during the process
	//	of selecting reads for normalisation
	if ((landscapeFilename) && !geneCounts.outputLandscape(landscapeFilename))
		exitFail("Unable to output landscape data to :", landscapeFilename);

	if (normalise)
	{
		coreParameterEstimation();

		elapsedTime("Parameter estimation complete");

		if (normaliseResultsFilename)
		{
			if (theModel == noModel)
			{
				bestModel = getBestModel();
				progMessage("Best model is ", bestModel);
				theModel = bestModel;
			}
			else
				progMessage("Model selected by command line is ", theModel);
		}
		else
		{
			progMessage("Model used is ", theModel);
		}

		getBias(theModel, bestResults[theModel].params[logValue], geneCounts.lengths[0], geneCounts.bias);

		if (normaliseResultsFilename)
		{
			//	Output the results of the mcmc analysis
			printResults();

			//	And then the bias predicted by all 6 models
			printBias();

			//	And then the counts and the bias for the genes themselves
			string filename = normaliseResultsFilename.replaceSuffix("_expression.txt");
			if (!geneCounts.outputGeneCounts(filename, 2, conv(theModel)))
				exitFail("Unable to output counts to :", filename);
		}
	}

	if(!geneCounts.outputGeneCounts(countsFilename,normalise?1:0,conv(theModel)))
		exitFail("Unable to output counts to :",countsFilename);

	elapsedTime("All results output");

	if (pauseAtEnd)
	{
		string test;
		cin >> test;
	}

	return EXIT_SUCCESS;
}



//	overlapCounts and chromosomeInfo would have been declared inside addRead as they are local to addRead.  
//	However this produces a C++ warning that the decorated name is too long, which can cause problems with debugging
struct overlapCounts
{
	size_t partial; // A count of the number of partial matches for this gene/region type combination
	size_t strict;	// A count of the number of strict matches for this gene/region type combination
	size_t length;  // The total length of the matches
	rna_pos_type RNAstartPos, RNAendPos;
	overlapCounts(size_t partial=0,size_t strict=0,size_t length = 0):partial(partial),strict(strict),length(length),RNAstartPos(99999999),RNAendPos(0){};
};

//	Contains information on the overlaps between a region of the read and featureRegions
//	The nested map is indexed first by gene /* and then by region type (e.g. exon) */ allowing data to be 
//	accumlated for multiple region types if required.
struct chromosomeGeneInfo: public map<string,overlapCounts> 
{
	size_t noMatch;
	size_t nSegments;
	chromosomeGeneInfo():noMatch(0),nSegments(0){};
};


void LiBiCount::addRead(const regionLists & segments,const featureFileEx & gtfData)
{
	if (segments.NH > 1)
	{
		geneCounts.count(notUnique)++;
		return;
	}
	else if (segments.qual < minqual)
	{
		geneCounts.count(lowQualString)++;
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
		//	For teh segments on each of the chromosomes (normally only one chromosome) get the featureRegions for the chromosome
		genomeFeatureRegions::const_iterator thisChromGtfRegions = gtfData.genomeGtfData.find(references[chromSegments.first].RefName);

		if (thisChromGtfRegions != gtfData.genomeGtfData.end())
		{
			chromosome = &references[chromSegments.first].RefName;

			if (outputFile.is_open() && location.empty())
			{
				location = stringEx(*chromosome,":",chromSegments.second.data.begin()->first);
			}

			//	Get the map of neds of featureRegions associated with the chromosome
			const chromosomeEndIndexMap & thisChromEndMap = gtfData.genomeEndIndex.at(references[chromSegments.first].RefName);

			if(thisChromEndMap.size())
			{

				//	And find the first one that finishes at or beyond the start of the first segment
				chromosomeEndIndexMap::const_iterator indirectIteratorStart = thisChromEndMap.lower_bound(chromSegments.second.data.begin()->second.start);

				//	And move back one to ensure we have the region that covers segment
	//			if (indirectIteratorStart != thisChromEndMap.begin())
				if (indirectIteratorStart == thisChromEndMap.end())
					indirectIteratorStart--;

				chromosomeFeatureData::iterator regionIterator = indirectIteratorStart->second;

				//	If this region overlaps any other regions then go to the one that starts the earliest.
				//	If there were no overlaps then default is the overlaps points to self
				regionIterator = *regionIterator->second.overlaps;


				//And now go through each of the segments
				for (auto & segment : chromSegments.second.data)
				{

					genes.nSegments++;
					vector<featureOverlap> overlaps;
					//	Trying out each of the featureRegions in turn to see if there is an overlap. 
					//	If there is then it gets added to the list of overlaps
					for (auto j = regionIterator; (j != thisChromGtfRegions->second.end()) && (j->first <= segment.second.end); j++)
					{
						//	Check for strand match
						if (!useStrand || ((j->second.strand == segment.second.strand) != reverseStrand))
						{
							j->second.checkOverlap(segment.second, overlaps);
						}

						regionIterator = *j->second.overlaps;
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
							string name, type;
							gtfId() {};
							gtfId(const featureOverlap & region) :name(region.geneName), type(region.featType) {};

						};

						class segRegion
						{
						public:
							size_t min, max;
							vector<gtfId> geneSet;
							segRegion() :min(INT_MAX), max(0) {};
							segRegion(const featureOverlap & region) :min(region.start), max(region.finish)
							{
								geneSet.emplace_back(region);
							};

						};

						//	
						featureOverlap & overlap = overlaps[0];

						_DBG(dbgFound = (overlap.geneName == "NM_001115075.1"));

						//	This will create an entry for the combination if it does not exist before, which is needed later on
						overlapCounts & GeneAttributeCombo1 = genes[overlap.geneName];
						GeneAttributeCombo1.RNAstartPos = min(overlap.RNAstartPos, GeneAttributeCombo1.RNAstartPos);
						GeneAttributeCombo1.RNAendPos = max(overlap.RNAendPos, GeneAttributeCombo1.RNAendPos);

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
							for (size_t i = 1; i < overlaps.size(); i++)
							{
								featureOverlap & overlap = overlaps[i];

								//	Create dummy entry for every region type that the read overlapped;

								overlapCounts & GeneAttributeCombo2 = genes[overlap.geneName];

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
										GeneAttributeCombo2.RNAstartPos = min(GeneAttributeCombo2.RNAstartPos, overlap.RNAstartPos);
										GeneAttributeCombo2.RNAendPos = max(GeneAttributeCombo2.RNAendPos, overlap.RNAendPos);
										newRegion = false;
									}
									else
									{
										newRegion = true;
										GeneAttributeCombo2.RNAstartPos = min(GeneAttributeCombo2.RNAstartPos, overlap.RNAstartPos);
										GeneAttributeCombo2.RNAendPos = max(GeneAttributeCombo2.RNAendPos, overlap.RNAendPos);

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
								for (size_t i = 1; i < segRegions.size(); i++)
								{
									bool found = false;
									for (auto & j : segRegions.at(i).geneSet)
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
									//	Use at() rather than [] as it is more efficient: assumes the entry is already in place
									genes.at(g.name).partial++;
									genes.at(g.name).length += size;
								}
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
	const string * mode = &blankString;		//The mode that selected the region

	rna_pos_type RNAstartPos = 0, RNAendPos = 0;

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
				if (gene.second.strict == genes.nSegments)
				{
					if (result)
					{
						if (countMode == intersect_all)
						{
							//	In intersect all we include all of the options, ie there will be multiple counts associated with
							//	a single fragment
							if (outputFile.is_open())
								outputFile.printEnd(*mode, *result, location, segments.name);
							geneCounts.count(*result)++;
							geneCounts.readPositionData[*result].positions[0].emplace_back(RNAstartPos);
							result = &gene.first;
							RNAstartPos = gene.second.RNAstartPos;
							RNAendPos = gene.second.RNAendPos;
						}
						else
						{
							//	If we have two strict matches then the result is ambigous, no need to look any further
							result = &ambiguousString;
							break;
						}
					}
					else
					{
						//	A strict match, keep looking as there may be more
						result = &gene.first;
						RNAstartPos = gene.second.RNAstartPos;
						RNAendPos = gene.second.RNAendPos;
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
			//	If the mode is intersection_nonempty we purposely move through to the
			//  intersection union case to see whether the read is a 'union' type of match
			//	so, no break between case statements
		}
	case intersect_union:
		{
			size_t Ngenes = genes.size();

			if (Ngenes == 1)
			{
				//	If we only match to one gene then the answer is simple
				result = &genes.begin()->first;
				RNAstartPos = genes.begin()->second.RNAstartPos;
				RNAendPos = genes.begin()->second.RNAendPos;
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
					//	length of the match within any of the segments is longest
					//	If multiple types are selected (exon/transcript) 
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
								//	Two genes with the same match length
								result = &ambiguousString;
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
		outputFile.printEnd(*mode,*result,location,segments.name);

	geneCounts.count(*result)++;

	rna_pos_type geneLen = geneCounts.length(*result);
	
	if (geneLen != 0)
	{
		if (segments.strands.size() == 1)
		{
			if (genomeDef.genes[*result].strand == '+')
			{
				if (segments.strands[0] == '+')
					geneCounts.readPositionData.at(*result).positions[0].emplace_back(max(min(RNAstartPos, geneLen - 1),(rna_pos_type)1));
				else
					geneCounts.readPositionData.at(*result).positions[1].emplace_back(max(min(RNAendPos,geneLen-1), (rna_pos_type)1));
			}
			else
			{
				if (segments.strands[0] == '+')
					geneCounts.readPositionData.at(*result).positions[1].emplace_back(max(geneLen - RNAstartPos + 1,(rna_pos_type)1));
				else
					geneCounts.readPositionData.at(*result).positions[0].emplace_back(max(geneLen - RNAendPos + 1,(rna_pos_type)1));
			}
		}
		else
		{
			static int messageCount = 0;
			// paired end, an entry for each end
			if (genomeDef.genes[*result].strand == '+')
			{
				if (segments.strands[0] == segments.strands[1])
				{
					geneCounts.readPositionData.at(*result).positions[0].emplace_back(max(min(RNAstartPos, geneLen - 1), (rna_pos_type)1));
					geneCounts.readPositionData.at(*result).positions[1].emplace_back(max(min(RNAendPos, geneLen - 1), (rna_pos_type)1));
				}
				else
				{
					if (messageCount++ < 30)
						optMessage("Mismatched paired end: ", segments.name, " Position: ",
							references[segments.data.begin()->first].RefName, ":", segments.data.begin()->second.data.begin()->first);
				}
			}
			else
			{
				if (segments.strands[0] == segments.strands[1])
				{
					geneCounts.readPositionData.at(*result).positions[1].emplace_back(max(geneLen - RNAstartPos + 1,(rna_pos_type)1));
					geneCounts.readPositionData.at(*result).positions[0].emplace_back(max(geneLen - RNAendPos + 1,(rna_pos_type)1));
				}
				else
				{
					if (messageCount++ < 30)
						optMessage("Mismatched paired end: ", segments.name, " Position: ",
							references[segments.data.begin()->first].RefName, ":", segments.data.begin()->second.data.begin()->first);
				}
			}
		}
	}
}

//	Increments counter as reads are read from the file
void LiBiCount::incBamCounter(const BamAlignment * ba,int size)
{
	if ((++bamCounter % REP_LEN) == 0)
		{
			stringEx msg(bamCounter, " BAM alignment record pairs processed.");
			if (ba)
			{
				if(ba->RefID >= 0)
					msg += _s(" cache size = ",size,"  ",references[ba->RefID].RefName,":",ba ->Position);
				else
					msg += _s(" cache size = ",size,"  unmapped read");
			}
			else if (size != -1)
				msg += _s(" cache reads processed = ",size);
			progMessage(msg);
		}
}

bool LiBiCount::AReadIsMapped(const BamAlignment & ba)
{
	if (!ba.IsMapped())
	{
		if (ba.IsPaired())
		{
			if (!ba.IsMateMapped())
			{
				if (ba.IsFirstMate())
					geneCounts.count(notAlignedString)++;
				return false;
			}
		}
		else
		{
			geneCounts.count(notAlignedString)++;
			return false;
		}
	}
	return true;
}

#define READ_BUFFER_SIZE 100

bool LiBiCount::processNameOrderedBamData()
{
	BamAlignment ba[READ_BUFFER_SIZE];
	bool used[READ_BUFFER_SIZE];

	bamCounter = 0;
	set<string> previousNames;

	bool OK = reader.GetNextAlignment(ba[0],false);

	while (OK)
	{
		string & name = ba[0].Name;
		_DBG(dbgFound = (name == BAMNAME);)

		int Nreads = 0;

		do {
			used[Nreads] = false;
			OK = reader.GetNextAlignment(ba[++Nreads],false);
		}
		while ((ba[Nreads].Name == name) && (OK) && (Nreads < READ_BUFFER_SIZE));

		//	This is a fairly direct implementation of the htseq logic in __init_.py line 570 onwards
		//	There are potential issues with the ability of the code to cope with datasets where there are multiple
		//	alignments
		if (AReadIsMapped(ba[0]))
		{
			if ((Nreads == 1) && (!ba[0].IsPaired()))
			{
				//
				addRead(regionLists(readData(move(ba[0])),name),genomeDef);
				incBamCounter();
			}
			else
			{
				for (int i = 0;i < Nreads;i++)
				{
					if (!used[i])
					{
						regionLists regions(readData(move(ba[i])),name);
						for (int j = i+1;j < Nreads;j++)
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

	//	A class for caching the reads for which we have not yet found a pair.  When it exceeds a certain size
	//	it gets saved to disk
	class readCacheClass : public map<string,vector<readData> >
	{
	public:
		void save(const stringEx & filename)
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
		_DBG(dbgFound = (ba.Name == BAMNAME);)

		//	The NH handling is complex is that there may be one NH (with NH = 1) at one end, and multiple NHs (with NH > 1) at the other
		//	The NH > 1 samples have to be used to either pair with the other, or to remove the NH = 1 sample 

		if (AReadIsMapped(ba))
		{
			if (ba.IsPaired())
			{
				//	Store reads in a cache so they can be paired up.
				stringEx mateIndex(ba.Name, "_",
#ifdef MATCH_USING_POSITION
					ba.IsMateMapped() ? stringEx(ba.MateRefID, ba.MatePosition) : "0", "_",
					ba.IsMapped() ? stringEx(ba.RefID, ba.Position) : "0", "_",
#endif
					(ba.IsMapped() && ba.IsMateMapped()) ? insertConv(-ba.InsertSize) : 0, "_",
					ba.IsFirstMate() ? "S" : "F");
				readCacheClass::iterator i = readCache.find(mateIndex);
				if (i == readCache.end())
				{
					stringEx thisIndex(ba.Name, "_",
#ifdef MATCH_USING_POSITION
						ba.IsMapped() ? stringEx(ba.RefID, ba.Position) : "0", "_",
						ba.IsMateMapped() ? stringEx(ba.MateRefID, ba.MatePosition) : "0", "_",
#endif
						(ba.IsMapped() && ba.IsMateMapped()) ? insertConv(ba.InsertSize) : 0, "_",
						ba.IsFirstMate() ? "F" : "S");

					i = readCache.find(thisIndex);
					if (i == readCache.end())
					{
						auto i = readCache.emplace(thisIndex, vector<readData>());
						i.first->second.emplace_back(move(ba));
					}
					else
						i->second.emplace_back(move(ba));
				}
				else
				{
					regionLists regions(i->second.front(), ba.Name);

					regions.combine(move(ba));

					addRead(regions, genomeDef);

					if (i->second.size() > 1)
						i->second.erase(i->second.begin());
					else
						readCache.erase(i);

					incBamCounter(&ba, readCache.size());

				}
			}
			else
			{
				regionLists regions(ba, ba.Name);
				addRead(regions, genomeDef);
			}
		}
		incBamCounter(&ba,readCache.size());

		if (readCache.size() > maxCacheSize)
		{
			optMessage("Outputting cache data ",cacheCounter+1);
			readCache.save({tempDirectory,"file",cacheCounter++});
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
				dbgFound = (tags[0] == BAMNAME);
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

//	Go through the reads in the file caches.  The read pairs associated with the same fragment will be in different files
//	and the files are ordered by read name so we only need to look at the next reads in each file to spot pairs that can be processed
//	as a pair
void LiBiCount::processCachedReads(size_t cacheFileCount)
{

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
				This[cacheReads[index].name][index].emplace_back(cacheReads[index].currentRead);

		}
	} reads;

	for (size_t i = 0;i < cacheFileCount;i++)
	{
		cacheReads[i].open(stringEx(tempDirectory,"file",i));
		reads[cacheReads[i].name][i].emplace_back(cacheReads[i].currentRead);
		for (int j = 0;(j < 10) && (cacheReads[i].readNext());j++)
		{
			reads[cacheReads[i].name][i].emplace_back(cacheReads[i].currentRead);
		}
	}

	int cacheReadCounter = 0;
	while (reads.size())
	{
		_DBG( dbgFound = (reads.begin()->first == BAMNAME);)

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

	vectorEx<size_t> cols{ {1,2,3,4,5}};


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

/*
	The methods associated with a cache file.  The cache is used for storing reads in 
	position ordered bam files until the mates turn up. Entries are added and removed 
	from the in memory cache until it reaches a certain size, when it is then written
	to disk
*/

LiBiCount::cacheEntry::~cacheEntry()
{
	close();
};


bool LiBiCount::cacheEntry::open(const std::string filename)
{
	file = new std::ifstream();
	file ->open(filename);
	if (!file ->is_open()) return false;
	fname = filename;
	readNext();
	return true;
}
bool LiBiCount::cacheEntry::readNext()
{
	if (file->eof())
		return false;
	currentRead.cigar.clear();
	std::string line;
	getline(*file,line);
	parseTsv(line,name, currentRead.refId, currentRead.position, currentRead.strand, 
		currentRead.cigar, currentRead.NH, currentRead.qual);
	return true;
}
void LiBiCount::cacheEntry::close()
{
	if (file)
	{
		file->close();
		delete (file);
		file = 0;
		remove(fname.c_str());
	}
}
