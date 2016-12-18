#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"
#include "mcmc.h"
#include "LogLiklihoods.h"

class LiBiNorm
{
public:
	LiBiNorm(size_t minModel,size_t maxModel,size_t Nruns) :
		minModel(minModel),maxModel(maxModel),Nruns(Nruns),
		outputMonteCarlo(false), outputConsolidated(false), outputNormalisation(false), singleModel(false) {};

	void mcmcThread(paramSet params, optionsType options, modelType model);
	int main(int argc, char **argv);
	bool core(size_t Nthreads,size_t Nsimu, int Ngenes);
	void printResults(const stringEx & outputFileName,const string & lastGene);
	void printNormalisation(const stringEx & outputFileName, const dataVec & lengths);

	transcriptDataMap transData;
	map<size_t, bestResult> bestResults;

private:
	size_t minModel,maxModel,Nruns;

	stringEx consFileName,outputFileName;

	dataType consData;
	bool outputMonteCarlo, outputConsolidated, outputNormalisation,singleModel;

	//	Four containers for the results.  The vector holds the results for each of the
	//	models.  For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished

	//	And these are for the full results within each chain
	vector<map<size_t,vector <dataVec > > >fullResultChain;
	vector<map<size_t,dataVec> >fullResultSSChain;

	dataVec RejectionRate;

	map<size_t, vector<string> > headers;

	//	Counts of the number of loops of each model
	map<size_t,int> threadLoopCounts;

	//	The results data
	map<size_t, multimap <double, dataVec *> > allOrderedResults;

};

#endif


