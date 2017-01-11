#ifndef LIBIOPTIMISER_H
#define LIBIOPTIMISER_H
#include "nelderMeadOptimiser.h"
#include "modelData.h"
#include "mcmc.h"

class LiBiOptimiser : public nelderMeadOptimiser
{
	optiDataType allOptiData;
	paramSet params;

public:

	optionsType * opts;
	const mcmcGeneData & geneData;
	LiBiOptimiser(const mcmcGeneData & geneData) : geneData(geneData){};

	dataVec getParams(modelType m,optionsType & options);


	virtual ErrorPair ErrorFunc();
	virtual void SaveResults(bool toFile) {};

};

#endif

