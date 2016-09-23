#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"

class LiBiNorm
{
	stringEx consFileName;
	transcriptDataMap transData;

	dataType consData;


	vector<vector <dataVec > > Chain;
	vector<dataVec> SSChain; 
	dataVec RejectionRate;

	map<size_t, vector<string> > headers;

	map<size_t,size_t> threadLoopCounts;

public:
	void mcmcThread(paramSet params, optionsType options, modelType model);
	int main(int argc, char **argv);
	int loadData();


};

#endif


