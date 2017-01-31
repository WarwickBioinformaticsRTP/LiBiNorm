#ifndef LIBINORM_H
#define LIBINORM_H

#include "LiBiNormCore.h"

class LiBiNorm : protected LiBiNormCore
{
public:
	LiBiNorm() : bestModel(noModel), outputFull(false) {};

	//	Used to identify param and param_dev values within bestResults
	enum
	{
		logValue = 0,
		absValue = 1,
		minLog = 2,
		maxLog = 3
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
	modelType bestModel;
private:
	bool outputFull;
	size_t nelderMeadCounter;

	//	These vector holds the full results for each of the models, which are needed for identifying
	//	the optimal parameters and the variation that is seen.
	//	For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished
	std::map<modelType,std::map<mcmcRunId, std::vector <dataVec > > >fullResultChain;
	std::map<modelType,std::map<mcmcRunId,dataVec> >fullResultSSChain;

	//	Counts of the number of mcmc runs that will be done for each model
	struct loop_counts { mcmcRunId requested, counter; };
	std::map<modelType, loop_counts> threadLoopCounts;

	//	Holds the results from multiple chains combined into a single ordered list
	std::map<modelType, std::multimap <double, dataVec *> > allOrderedResults;

};

#endif


