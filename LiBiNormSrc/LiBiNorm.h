#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"

class LiBiNorm
{
	stringEx consFileName,outputFileName;
	transcriptDataMap transData;

	dataType consData;
	bool fullOutputMode,singleModel;


	vector<map <size_t,dataVec > > Chain;
	vector<map <size_t, VEC_DATA_TYPE> > SSChain;

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


