#ifndef LIBINORM_H
#define LIBINORM_H

#include <map>
#include "stringEx.h"
#include "containerEx.h"
#include "GeneCountData.h"
#include "ModelData.h"


#define DEF_MAX_READS_FOR_PARAM_ESTIMATION 1000000
//#define MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 20000
#define MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 200000000
#define MAX_READS_GENE 100
#define DEF_THREADS 3
#define MCMC_ITERATIONS 2000
#define DEFAULT_MODEL ModelBD
#define NUMBER_OF_MCMC_RUNS 3
#define DEFAULT_NORMALISATION_GENE_LENGTH 1000
#define MAX_GENE_LENGTH_FOR_NORM_PLOT 20000
#define DEFAULT_COUNT_MODE intersect_union
#define DEFAULT_FEATURE_TYPE_EXON "exon" 
#define DEFAULT_GTF_ID_ATTRIBUTE "gene_id"
#define DEFAULT_GFF_ID_ATTRIBUTE "Genbank"

//	Use this to add the mode which creates fastq files based on bam files with artificial problems
// #define MAKE_FASTQ_MODE

//	Use this to add the mode where duplicates in bam files can be removed
// #define DEDUP_MODE

//	Use this mode to run a model with specific parameters.  This affects how the values are set
//	in ModelParameters.cpp and also ensures that the first N reads are use in data are loaded
//	in void rnaPosVec::selectAtMost(size_t s) in GeneCountData.cpp
//#define PRESET_VALUES {-0.060374433,	1.944021977,	-4.434270212,	-3.367865765,	0.877799}

//	For identifying the mcmc run for a specific model.  Numbered from 1.
typedef size_t mcmcRunId;

class LiBiNormCore
{
protected:
	LiBiNormCore() :normalise(true), pauseAtEnd(false),
		theModel(noModel),
		maxReads(DEF_MAX_READS_FOR_PARAM_ESTIMATION),
		Nthreads(DEF_THREADS),
		Nsimu(MCMC_ITERATIONS),
		Nruns(NUMBER_OF_MCMC_RUNS),
		NrunsOtherModels(0)
	{};

	void helpCommon();
	bool commandParseCommon(int & ni, int argc, char **argv);

	bool normalise, pauseAtEnd;
	modelType theModel;
	size_t maxReads, Nthreads, Nsimu;
	mcmcRunId Nruns,NrunsOtherModels;

	stringEx landscapeFilename, normaliseResultsFilename, countsFilename;
};

typedef std::map<modelType, std::vector<std::string> > headerType;

class LiBiNorm : protected LiBiNormCore
{
public:
	LiBiNorm() : bestModel(noModel),outputFull(false) 
	{
		headers = getHeaders();
	};

	int main(int argc, char **argv);
	void mcmcThread(optionsType options);
protected:
	bool coreParameterEstimation();
	modelType getBestModel();

	//	For outputting results
	void printResults();
	void printBias();
	void printAllMcmcRunData();
	void printConsolidatedMcmcRunData();

protected:
	std::map<modelType, bestResult> bestResults;
	GeneCountData geneCounts;
	modelType bestModel;
private:
	bool outputFull;

	//	The specific data that will be used for the mcmc parameter determination
	mcmcGeneData geneData;

	//	These vector holds the full results for each of the models, which are needed for identifying
	//	the optimal parameters and the variation that is seen.
	//	For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished
	std::map<modelType,std::map<mcmcRunId, std::vector <dataVec > > >fullResultChain;
	std::map<modelType,std::map<mcmcRunId,dataVec> >fullResultSSChain;

	headerType headers;

	//	Counts of the number of mcmc runs that will be done for each model
	struct loop_counts { mcmcRunId requested, counter; };
	std::map<modelType, loop_counts> threadLoopCounts;

	//	Holds the results from multiple chains combined into a single ordered list
	std::map<modelType, std::multimap <double, dataVec *> > allOrderedResults;

};

#endif


