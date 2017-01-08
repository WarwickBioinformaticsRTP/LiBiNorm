#ifndef LIBINORM_H
#define LIBINORM_H

#include <map>
#include "stringEx.h"
#include "containerEx.h"
#include "GeneCountData.h"
#include "Options.h"
#include "ModelData.h"

//	The original MATLAB code had an error in setting the initial values for mcmc runs which this 
//	simulates (ModelData.cpp)
// #define SIMULATE_MATLAB_BUG

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


