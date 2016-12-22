#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"
#include "mcmc.h"
#include "LogLiklihoods.h"


#define DEF_MAX_READS_FOR_PARAM_ESTIMATION 1000000
#define MAX_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 20000
#define MAX_READS_GENE 100
#define DEF_THREADS 3
#define MCMC_ITERATIONS 2000
#define N_MODELS 6
#define DEFAULT_MODEL 6
#define NUMBER_OF_MCMC_RUNS 3
#define DEFAULT_NORMALISATION_GENE_LENGTH 1000
#define MAX_GENE_LENGTH_FOR_NORM_PLOT 20000
#define DEFAULT_COUNT_MODE intersect_union
#define DEFAULT_FEATURE_TYPE_EXON "exon" 
#define DEFAULT_GTF_ID_ATTRIBUTE "gene_id"
#define DEFAULT_GFF_ID_ATTRIBUTE "Genbank"


class LiBiNormCore
{
protected:
	LiBiNormCore() :normalise(false), 
		maxReads(DEF_MAX_READS_FOR_PARAM_ESTIMATION),
		Nthreads(DEF_THREADS),
		Nsimu(MCMC_ITERATIONS),
		Nruns(NUMBER_OF_MCMC_RUNS),
		NrunsOtherModels(0)
	{};

	void helpCommon();
	bool commandParseCommon(int & ni, char **argv);

	bool normalise;
	size_t maxReads, Nthreads,Nsimu,Nruns,NrunsOtherModels;

	stringEx landscapeFilename, normaliseResultsFilename, countsFilename;
};


class LiBiNorm : protected LiBiNormCore
{
public:
	LiBiNorm() :
		theModel(DEFAULT_MODEL), outputFull(false) {};

	void mcmcThread(paramSet params, optionsType options, modelType model);
	int main(int argc, char **argv);
	bool coreParameterEstimation();
	void printResults(const string & lastGene);
	void printBias();
	void printAllMcmcRunData();
	void printConsolidatedMcmcRunData();

//	transcriptDataMap transData;
	map<size_t, bestResult> bestResults;

protected:
	GeneCountData geneCounts;

private:
	size_t theModel;
	bool outputFull;

	dataType consData;

	//	These vector holds the full results for each of the models, which are needed for identifying
	//	the optimal parameters and the variation that is seen.
	//	For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished
	vector<map<size_t,vector <dataVec > > >fullResultChain;
	vector<map<size_t,dataVec> >fullResultSSChain;

	dataVec RejectionRate;

	map<size_t, vector<string> > headers;

	//	Counts of the number of mcmc runs that will be done for each model
	map<size_t,int> threadLoopCounts;

	//	The results data
	map<size_t, multimap <double, dataVec *> > allOrderedResults;

};

#endif


