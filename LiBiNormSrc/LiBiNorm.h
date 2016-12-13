#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"

extern bool verbose;

#define STORE_ENDPOINTS

class LiBiNorm
{
	stringEx consFileName,outputFileName;
	transcriptDataMap transData;

	dataType consData;
	bool fullOutputMode,singleModel;

	//	Four containers for the results.  The vector holds the results for each of the
	//	models.  For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished

	//	These two are for the endpoints of the chain
#ifdef STORE_ENDPOINTS
	vector<map <size_t,dataVec > > Chain;
	vector<map <size_t, VEC_DATA_TYPE> > SSChain;
#endif
	//	And these are for the full results within each chain
	vector<map<size_t,vector <dataVec > > >fullResultChain;
	vector<map<size_t,dataVec> >fullResultSSChain;

	dataVec RejectionRate;

	map<size_t, vector<string> > headers;

	//	Counts of the number of loops of each model
	map<size_t,int> threadLoopCounts;

public:
	void mcmcThread(paramSet params, optionsType options, modelType model);
	int main(int argc, char **argv);
	int loadData();


};

#endif


