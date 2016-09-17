#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"

class LiBiNorm
{
	stringEx consFileName;
	transcriptDataMap transData;

	dataType consData;


	vector<dataVec > Chain;
	vector<dataVec> SSChain; 
	dataVec RejectionRate;

	size_t threadLoopCount;

public:
	void mcmcThread(paramSet params, optionsType options, modelType model);
	int main(int argc, char **argv);
	int loadData();


};

#endif


